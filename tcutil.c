/* tcutil.c:  Wrapper to hide all the termcap vilesomeness */

#ifdef _WIN32

#include <stdio.h>
#include <stdlib.h>

#include <windows.h>

#include "tcutil.h"

HANDLE h_stdout;

int tcutil_init (void) {
  CONSOLE_SCREEN_BUFFER_INFO inf;

  /* Ensure stdio buffers are flushed */
  fflush(stdout);

  /* Get handle */
  h_stdout = GetStdHandle(STD_OUTPUT_HANDLE);
  if(h_stdout == INVALID_HANDLE_VALUE)
    return(0);

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

void tcutil_back (void) {
  CONSOLE_SCREEN_BUFFER_INFO inf;

  if(GetConsoleScreenBufferInfo(h_stdout, &inf)) {
    if(inf.dwCursorPosition.X > 0) {
      inf.dwCursorPosition.X--;
      SetConsoleCursorPosition(h_stdout, inf.dwCursorPosition);
    }
  }
}

void tcutil_up (void) {
  CONSOLE_SCREEN_BUFFER_INFO inf;

  if(GetConsoleScreenBufferInfo(h_stdout, &inf)) {
    if(inf.dwCursorPosition.Y > 0) {
      inf.dwCursorPosition.Y--;
      SetConsoleCursorPosition(h_stdout, inf.dwCursorPosition);
    }
  }
}

#define FG_GREY (FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE)
#define BG_GREY (BACKGROUND_RED | BACKGROUND_GREEN | BACKGROUND_BLUE)

void tcutil_attr (int a) {
  CONSOLE_SCREEN_BUFFER_INFO inf;

  if(GetConsoleScreenBufferInfo(h_stdout, &inf)) {
    switch(a) {
      /* These are mutually exclusive, so we just set it */
    case ATTR_NORM:
    case ATTR_BLINK:
    case ATTR_HALF:
      inf.wAttributes = FG_GREY;
      break;
    case ATTR_BOLD:
      inf.wAttributes = FG_GREY | FOREGROUND_INTENSITY;
      break;
    case ATTR_REVERSE:
    case ATTR_STANDOUT:  /* XXX - this is not really what ANSI does */
      inf.wAttributes = BG_GREY;
      break;
    case ATTR_UNDERLINE:
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

  /* Look up backspace */
  p = tgetstr("bc", &bp);
  if(p)
    BC = p;
  else {
    /* Check for backspace capability */
    if(tgetflag("bs"))
      BC = (char *) NULL;
    else  /* not enough cursor motion capability */
      return(0);
  }

  /* Look up "move up a line" */
  p = tgetstr("up", &bp);
  if(p)
    UP = p;
  else
    return(0);

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

void tcutil_back (void) {
  if(BC)
    tputs(BC, 1, putchar);
  else
    putchar('\b');
}

void tcutil_up (void) {
  tputs(UP, 1, putchar);
}

void tcutil_attr (int a) {
  if(attr[a])
    tputs(attr[a], 1, putchar);
}

#endif  /* _WIN32 */
