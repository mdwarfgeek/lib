/* tcutil.c:  Primitive terminal manipulation functions. */

#ifdef _WIN32

#include <stdio.h>
#include <stdlib.h>

#include <windows.h>

#include "tcutil.h"

static HANDLE h_stdout;
static unsigned char reversed = 0;

int tcutil_init (void) {
  CONSOLE_SCREEN_BUFFER_INFO inf;

  /* Ensure stdio buffers are flushed */
  fflush(stdout);

  /* Get handle */
  h_stdout = GetStdHandle(STD_OUTPUT_HANDLE);
  if(h_stdout == INVALID_HANDLE_VALUE)
    return(0);

  /* Init reverse video flag */
  reversed = 0;

  /* Get info */
  if(!GetConsoleScreenBufferInfo(h_stdout, &inf))
    return(0);
  
  /* Turn off buffering in stdio */
  setvbuf(stdout, (char *) NULL, _IONBF, 0);

  return(1);
}

void tcutil_winsize (int *ncols, int *nrows) {
  CONSOLE_SCREEN_BUFFER_INFO inf;

  if(GetConsoleScreenBufferInfo(h_stdout, &inf)) {
    *ncols = inf.dwSize.X;
    *nrows = inf.dwSize.Y;
  }
  else {
    *ncols = 0;
    *nrows = 0;
  }
}

void tcutil_goto (int x, int y) {
  COORD cp;

  cp.X = x;
  cp.Y = y;

  SetConsoleCursorPosition(h_stdout, cp);
}

void tcutil_back (int n) {
  CONSOLE_SCREEN_BUFFER_INFO inf;

  if(GetConsoleScreenBufferInfo(h_stdout, &inf)) {
    if(inf.dwCursorPosition.X > 0) {
      if(inf.dwCursorPosition.X > n)
        inf.dwCursorPosition.X -= n;
      else
        inf.dwCursorPosition.X = 0;

      SetConsoleCursorPosition(h_stdout, inf.dwCursorPosition);
    }
  }
}

void tcutil_up (int n) {
  CONSOLE_SCREEN_BUFFER_INFO inf;

  if(GetConsoleScreenBufferInfo(h_stdout, &inf)) {
    if(inf.dwCursorPosition.Y > 0) {
      if(inf.dwCursorPosition.Y > n)
        inf.dwCursorPosition.Y -= n;
      else
        inf.dwCursorPosition.Y = 0;

      SetConsoleCursorPosition(h_stdout, inf.dwCursorPosition);
    }
  }
}

#define FG_GREY (FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE)
#define BG_GREY (BACKGROUND_RED | BACKGROUND_GREEN | BACKGROUND_BLUE)

void tcutil_attr (int a) {
  CONSOLE_SCREEN_BUFFER_INFO inf;
  WORD tmp;

  if(GetConsoleScreenBufferInfo(h_stdout, &inf)) {
    switch(a) {
    case ATTR_NORM:
      /* Reset to "normal" */
      inf.wAttributes = FG_GREY;
      reversed = 0;
      break;
    case ATTR_BLINK:
    case ATTR_HALF:
      /* These do nothing, we don't have this capability. */
      break;
    case ATTR_BOLD:
      /* Add bold to what's there using intensity flag. */
      if(inf.wAttributes & FG_GREY)  /* any foreground, intensify that */
        inf.wAttributes |= FOREGROUND_INTENSITY;
      else  /* otherwise the background */
        inf.wAttributes |= BACKGROUND_INTENSITY;
      break;
    case ATTR_REVERSE:
    case ATTR_STANDOUT:
      /* Most terminals don't revert to normal video when sent
         multiple reverses.  We emulate this behaviour using a
         flag set at the first reverse to prevent more from
         being executed. */
      if(!reversed) {
        /* Swap foreground and background colours and flags.
           Done the hard way to avoid depending on the exact
           values and storage format for the flags. */
        tmp = 0;
        
        if(inf.wAttributes & FOREGROUND_RED)
          tmp |= BACKGROUND_RED;
        if(inf.wAttributes & FOREGROUND_GREEN)
          tmp |= BACKGROUND_GREEN;
        if(inf.wAttributes & FOREGROUND_BLUE)
          tmp |= BACKGROUND_BLUE;
        if(inf.wAttributes & FOREGROUND_INTENSITY)
          tmp |= BACKGROUND_INTENSITY;
        
        if(inf.wAttributes & BACKGROUND_RED)
          tmp |= FOREGROUND_RED;
        if(inf.wAttributes & BACKGROUND_GREEN)
          tmp |= FOREGROUND_GREEN;
        if(inf.wAttributes & BACKGROUND_BLUE)
          tmp |= FOREGROUND_BLUE;
        if(inf.wAttributes & BACKGROUND_INTENSITY)
          tmp |= FOREGROUND_INTENSITY;

        inf.wAttributes = tmp;
        reversed = 1;
      }
      break;
    case ATTR_UNDERLINE:
      /* Underline is cyan (or more generally, no red). */
      inf.wAttributes &= ~(FOREGROUND_RED | BACKGROUND_RED);
      break;
    }

    SetConsoleTextAttribute(h_stdout, inf.wAttributes);
  }
}

#else

#include <sys/types.h>
#include <sys/ioctl.h>
#include <termios.h>
#include <unistd.h>

#include <stdio.h>
#include <stdlib.h>
#include <curses.h>
#include <term.h>

#include "tcutil.h"

extern char PC;
extern char *BC;
extern char *UP;
extern short ospeed;

#define OP_BC    0
#define OP_UP    1
#define OP_BCN   2
#define OP_UPN   3
#define OP_CM    4
#define NUM_OPS  5

static char *ops_tc[NUM_OPS] = { "bc", "up", "BC", "UP", "cm" };
static char *ops[NUM_OPS];
static unsigned char can_bs;

static char *attr_tc[NUM_ATTR] = { "me", "mb", "md", "mh", "mr", "so", "us" };
static char *attr[NUM_ATTR];

static char termcapbuf[2048];
static char termentbuf[2048];

int tcutil_init (void) {
  char *termtype, *p, *bp;

  int a, rv;

  /* Check we have a real terminal */
  if(!isatty(0))
    return(0);

  /* Get terminal type */
  termtype = getenv("TERM");
  if(!termtype)
    return(0);

  /* Get termcap */
  rv = tgetent(termcapbuf, termtype);
  if(rv != 1)
    return(0);

  bp = termentbuf;

  /* Look up padding */
  p = tgetstr("pc", &bp);
  if(p)
    PC = *p;
  else
    PC = 0;

  /* Can the terminal backspace?  We may need this. */
  can_bs = tgetflag("bs");

  /* Set BC and UP to NULL so we can use tgoto to emit other caps
     needing arguments. */
  BC = (char *) NULL;
  UP = (char *) NULL;

  /* Look up operations */
  for(a = 0; a < NUM_OPS; a++) {
    p = tgetstr(ops_tc[a], &bp);
    if(p)
      ops[a] = p;
    else
      ops[a] = (char *) NULL;
  }

  /* Look up attributes */
  for(a = 0; a < NUM_ATTR; a++) {
    p = tgetstr(attr_tc[a], &bp);
    if(p)
      attr[a] = p;
    else
      attr[a] = (char *) NULL;
  }

  /* Turn off buffering */
  setvbuf(stdout, (char *) NULL, _IONBF, 0);

  /* Should really turn off output processing
   * and arrange for restoration on exit, but that
   * would mean we have to handle all the other
   * stuff (like newlines) ourselves too.
   */

  /* Success! */
  return(1);
}

void tcutil_winsize (int *ncols, int *nrows) {
  int rv;
  struct winsize wsz;

  /* Try the ioctl first */
  rv = ioctl(0, TIOCGWINSZ, &wsz);
  if(rv == 0) {
    *ncols = wsz.ws_col;
    *nrows = wsz.ws_row;
  }
  else {
    /* Try termcap */
    *ncols = tgetnum("co");
    *nrows = tgetnum("li");
  }
}

void tcutil_goto (int x, int y) {
  char *s;

  if(ops[OP_CM]) {
    /* We need to set these to use tgoto */
    BC = ops[OP_BC];
    UP = ops[OP_UP];

    /* Get the string */
    s = tgoto(ops[OP_CM], x, y);

    /* Check if we got anything, and emit */
    if(s)
      tputs(s, 1, putchar);

    /* Restore */
    BC = (char *) NULL;
    UP = (char *) NULL;
  }
}

void tcutil_back (int n) {
  int i;
  char *s;

  if(n > 1 && ops[OP_BCN]) {
    s = tgoto(ops[OP_BCN], 0, n);
    if(s)
      tputs(s, 1, putchar);
  }
  else if(ops[OP_BC]) {
    for(i = 0; i < n; i++)
      tputs(ops[OP_BC], 1, putchar);
  }
  else if(can_bs) {
    for(i = 0; i < n; i++)
      putchar('\b');
  }
}

void tcutil_up (int n) {
  int i;
  char *s;

  if(n > 1 && ops[OP_UPN]) {
    s = tgoto(ops[OP_UPN], 0, n);
    if(s)
      tputs(s, 1, putchar);
  }
  else if(ops[OP_UP]) {
    for(i = 0; i < n; i++)
      tputs(ops[OP_UP], 1, putchar);
  }
}

void tcutil_attr (int a) {
  if(attr[a])
    tputs(attr[a], 1, putchar);
}

#endif  /* _WIN32 */
