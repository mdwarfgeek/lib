#ifndef CNF_H
#define CNF_H

struct cnf_field {
  char *section;      /* [section] name */
  char *key;          /* keyword */
  enum {
    CNF_FIELD_INT,
    CNF_FIELD_FLOAT,
    CNF_FIELD_DOUBLE,
    CNF_FIELD_DBL60,
    CNF_FIELD_STRING,
    CNF_FIELD_PATH
  } type;             /* code for field type */
  void *ptr;          /* pointer to receive data */
  int len;            /* length of array (for string) */
  int found;          /* (output): seen? */
};

int cnf_read (char *filename,
	      struct cnf_field *fields, int nfields,
	      char *errstr);

#endif  /* CNF_H */
