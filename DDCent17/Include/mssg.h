/*****************************************************************************
**
**  FILE:    mssg.h
**
**  AUTHOR:  Melannie Hartman  9/7/93 - 8/21/96
**
**  Updates:
**    WFPS effect on N2O flux                  KLK 7 Aug 13
**    buffer output for internal redirection   KLK
**
*****************************************************************************/

#include <stdio.h>

char msgbfr[256];

void mssgfil(char filnam[256]);
void closmssgfil();
void mssg(char buffer[256]);
void strout(char buffer[256]);
void abortmssg(char buffer[256]);
/* abortmssg with translated errno; use only if errno is known to be valid */
void aborterr(char buffer[256]);

