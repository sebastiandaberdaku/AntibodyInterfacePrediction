/*
 * disclaimer.h
 *
 *  Created on: 29/ago/2014
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_UTILS_DISCLAIMER_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_UTILS_DISCLAIMER_H_

#include <iostream>

/**
 * Author
 */
#define AUTHOR "Sebastian Daberdaku"
/**
 * email
 */
#define EMAIL "sebastian.daberdaku@dei.unipd.it"
/**
 * 3-clause BSD License.
 */
#define DISCLAIMER \
	cout <<\
"Copyright (c) 2017, Sebastian Daberdaku - \n\
DIPARTIMENTO DI INGEGNERIA DELL'INFORMAZIONE - \n\
Università degli Studi di Padova.\n\
All rights reserved.\n\n\
\
Redistribution and use in source and binary forms, with or without\n\
modification, are permitted provided that the following conditions are met:\n\
  * Redistributions of source code must retain the above copyright\n\
    notice, this list of conditions and the following disclaimer.\n\
  * Redistributions in binary form must reproduce the above copyright\n\
    notice, this list of conditions and the following disclaimer in the\n\
    documentation and/or other materials provided with the distribution.\n\
  * The name of the author may not be used to endorse or promote products\n\
    derived from this software without specific prior written permission.\n\n\
\
THIS SOFTWARE IS PROVIDED BY THE AUTHOR ''AS IS'' AND ANY EXPRESS OR\n\
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES\n\
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.\n\
IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,\n\
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT\n\
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,\n\
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY\n\
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT\n\
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF\n\
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n\n\n";
/**
 * Name of the program
 */
#define PROGRAM_NAME "Paratope prediction - Test samples generator.\n"/**
 * Program version
 */
#define PROGRAM_VERSION "1.0"
/**
 * Program information.
 */
#define PROGRAM_INFO \
	cout <<\
"\n"<<\
"**************************************************\n"<<\
" Paratope prediction - Test samples generator.\n"<<\
" Interface prediction program for Ab-Ag\n"<<\
" complexes based on surface-patch descriptors.\n"<<\
" Program: "<<argv[0]<<" in source "<<__FILE__<<",\n"<<\
" v"<<PROGRAM_VERSION<<", compiled on "<<__DATE__<<" at "<<__TIME__<<".\n"<<\
" Author: "<<AUTHOR<<"\n"<<\
" Contact: "<<EMAIL"\n"<<\
"**************************************************\n\n";


#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_UTILS_DISCLAIMER_H_ */
