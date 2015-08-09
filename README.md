# cl-hmm
A small Hidden Markov Model package.

## Installation
With ASDF:  
If using ASDF 3.1.2 or later, you can move the source tree to ~/common-lisp/ and print:

    CL-USER> (asdf:load-system :cl-hmm)

or

    CL-USER> (require :cl-hmm)

if you are using recent version of ABCL, Clozure CL, CMUCL, ECL, CLISP, MKCL
and SBCL, which support this hook of ASDF.  
If using ASDF with version earlier than 3.1.2, move the source tree to ~/.local/share/common-lisp/source
and do the same thing.

Without ASDF:  
Just load the files into your REPL, then we're done.
