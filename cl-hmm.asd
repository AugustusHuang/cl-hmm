(in-package :asdf-user)
;(in-package :cl-user)

;(defpackage :cl-hmm-system
 ; (:use :cl :asdf))

;(in-package :cl-hmm-system)

(defsystem "cl-hmm"
  :description "cl-hmm: a small Hidden Markov Model system."
  :version "0.0.1"
  :author "Augustus Huang <augustushwang@gmail.com>"
  :components ((:file "cl-hmm")))
