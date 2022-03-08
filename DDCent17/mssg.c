#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h> /*Included for dup(2) and write(2)*/
#include "mssg.h"

static int oldstdout =-1;

FILE * termFile;

  /* char sq[]  = "set quiet";*/

void mssgfil(char filnam[256]) {
    fflush(stdout); /* write any STDOUT buffered to terminal*/
    oldstdout = dup(STDOUT_FILENO); /* save the existing STDOUT*/
    termFile = freopen(filnam, "w", stdout);    /* reroute STDOUT*/

    /* Error message and terminate if we can't redirect output */
    if (termFile == NULL) {
      printf  ("\nERROR: unable to open terminal %s\nAbnormal Termination\n",filnam);
    }
    /*else {printf  ("\nopen terminal %s\n",filnam);}*/

    return;
}

void closmssgfil() {
    fflush(stdout); /* write any STDOUT buffered to the file STDOUT*/
    close(STDOUT_FILENO); /*close the file */
    dup2(oldstdout, STDOUT_FILENO); /* copy the original terminal to STDOUT*/
    return;
}

void mssg(char buffer[256])
{
  int ret;
  ret = printf ("%s\n",buffer);
  /* if the stdout fails write the message to stderr and quit */
  if(ret < 0) {
    fprintf(stderr, "ERROR writing '%s' to terminal\nAbnormal Termination\n",buffer);
    exit(1);
  }
  return;
}

void strout(char buffer[256])
{
  int ret;
  ret = printf ("%s",buffer);
  /* if the stdout fails write the message to stderr and quit */
  if(ret < 0) {
    fprintf(stderr, "ERROR writing '%s' to terminal\nAbnormal Termination\n",buffer);
    exit(1);
  }
  return;
}

void aborterr(char buffer[256]) {
  int ret;

   if (errno>0) {ret = printf("ERROR: %s; %s\nAbnormal Termination\n", buffer, strerror(errno));}
   else         {ret = printf("ERROR: %s\nAbnormal Termination\n", buffer);}

  /*   ret = printf  ("%s\nAbnormal Termination\n",buffer);
if stdout fails try to stderr */
  if(ret < 0) {
   if (errno>0) {ret = fprintf(stderr,"stderr ERROR: %s; %s\nAbnormal Termination\n", buffer, strerror(errno));}
   else         {ret = fprintf(stderr,"stderr ERROR: %s\nAbnormal Termination\n", buffer);}
  }
  /* If the output is redirected; flush buffers and write "abnormal termination to STDERR*/
  /* if(oldstdout >= 0) {fflush((char *)NULL); write(oldstdout, "Abnormal Termination\n", 21);} */
  if(oldstdout >= -1) {fflush(stdout); fprintf(stderr, "Abnormal Termination\n");}
  exit(1);
}

void abortmssg(char buffer[256]) {
  int ret;
  ret = printf  ("ERROR: %s\nAbnormal Termination\n",buffer);
  /* if stdout fails try to stderr */
  if(ret < 0) {
    fprintf(stderr, "ERROR writing '%s' to terminal\nAbnormal Termination\n",buffer);
  }
  /* If the output is redirected; flush buffers and write "abnormal termination" to STDERR*/
  if(oldstdout >= -1) {fflush(stdout); fprintf(stderr, "Abnormal Termination\n");}
  exit(1);
}

/*
   strncmp	 Compare characters of two strings (function)
   memcmp	 Compare two blocks of memory (function)
   strrchr	 Locate last occurrence of character in string (function)
   strspn	 Get span of character set in string (function)
*/
