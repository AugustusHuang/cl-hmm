;;;; The MIT License (MIT)

;;;; Copyright (c) 2015 Augustus Huang

;;;; Permission is hereby granted, free of charge, to any person obtaining a copy
;;;; of this software and associated documentation files (the "Software"), to deal
;;;; in the Software without restriction, including without limitation the rights
;;;; to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
;;;; copies of the Software, and to permit persons to whom the Software is
;;;; furnished to do so, subject to the following conditions:

;;;; The above copyright notice and this permission notice shall be included in all
;;;; copies or substantial portions of the Software.

;;;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
;;;; IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
;;;; FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
;;;; AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
;;;; LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
;;;; OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
;;;; SOFTWARE.

;;;; Common Lisp Hidden Markov Model Algorithms Package
;;;; Date: May 22 2015

;;; Tested on SBCL and CLISP
(defpackage :cl-hmm-algorithms
  (:nicknames :hmm :cl-hmm)
  (:use :cl)
  (:export :get-states
	   :get-observations
	   :get-initial
	   :get-transition
	   :get-measure
	   :get-gamma
	   :get-xi
	   :init-hmm
	   :init-from-file
	   :print-to-file
	   :forward
	   :backward
	   :viterbi))

(in-package :cl-hmm-algorithms)

;;; Matrix type, in fact a 2-dimensional array.
(deftype matrix (&optional type x y)
  "Transition square matrix and measure rectangle matrix."
  `(array ,type (,x ,y)))

;;; In order to pretty print an array, implement it.
(defun pprint-array (stream array)
  "Pretty print routine to print a nested vector."
  (loop with first-time = t
     for x across array unless first-time
     do (write-char #\Space stream) end
     do (princ x stream)
       (setf first-time nil)))

;;; And so with a matrix...
(defun pprint-matrix (stream matrix)
  "Pretty print routine to print a matrix."
  (loop for i below (car (array-dimensions matrix)) do
       (loop for j below (cadr (array-dimensions matrix)) do
	    (let ((cell (aref matrix i j)))
	      (format stream "~a " cell)))
       (format stream "~%")))

;;; HMM state will be made up with:
;;; states --- internal states
;;; measures --- measurements or observations by external things
;;; initial --- initial distribution of internal states
;;; transition-matrix --- probabilities of transitions from a internal state to
;;; another
;;; measure-matrix --- probabilities of an internal state measured with output
;;; a specified observation
(defstruct hmm-state
  (states 0 :type fixnum)
  (measures 0 :type fixnum)
  (initial #(0) :type vector)
  (transition-matrix #2A((0)) :type matrix)
  (measure-matrix #2A((0)) :type matrix))

;;; Wrapper functions.
(declaim (inline get-states
		 get-observations
		 get-initial
		 get-transition
		 get-measure
		 (setf get-states)
		 (setf get-observations)
		 (setf get-initial)
		 (setf get-transition)
		 (setf get-measure)))

(defun get-states (hstate)
  (hmm-state-states hstate))

(defun (setf get-states) (val hstate)
  (setf (hmm-state-states hstate) val))

(defun get-observations (hstate)
  (hmm-state-measures hstate))

(defun (setf get-observations) (val hstate)
  (setf (hmm-state-measures hstate) val))

(defun get-initial (hstate)
  (hmm-state-initial hstate))

(defun (setf get-initial) (val hstate)
  (setf (hmm-state-initial hstate) val))

(defun get-transition (hstate)
  (hmm-state-transition-matrix hstate))

(defun (setf get-transition) (val hstate)
  (setf (hmm-state-transition-matrix hstate) val))

(defun get-measure-matrix (hstate)
  (hmm-state-measure-matrix hstate))

(defun (setf get-measure-matrix) (val hstate)
  (setf (hmm-state-measure-matrix hstate) val))

;;; Helpers, to generate a new HMM state from file.
(defun list-array (lst)
  "Wrapper function in order to build an array from a list."
  (make-array (length lst) :initial-contents lst))

(defun list-matrix (lst)
  "Wrapper function in order to build a matrix from a nested list."
  (make-array (list (length lst)
		    (length (first lst)))
	      :initial-contents lst))

;;; Only when we need to coerce ratio into double-float...
;;; Or I will keep floating-point ratio...
(defun parse-decimal (str &key (start 0) (end (length str)))
  "Read string with format 'abc.efg' and parse it."
  (let ((sum 0)
	(digit 0)
	(in-float-p nil)
	(float-depth 0))
    (do ((point start (1+ point)))
	((>= point end))
      (cond ((and (>= (char-code (char str point)) (char-code #\0))
		  (<= (char-code (char str point)) (char-code #\9)))
	     ;; We are meeting a number.
	     (progn
	       (setf digit (- (char-code (char str point)) 48)
		     sum (+ (* sum 10) digit))
	       ;; If we are in the floating part...
	       (if in-float-p
		   (incf float-depth))))
	    ((= (char-code (char str point)) (char-code #\.))
	     (if in-float-p
		 ;; Can't be in floating part twice.
		 (error "Invalid number format.")
		 (setf in-float-p t)))
	    (t
	     (error "Invalid number format."))))
    (/ sum (expt 10 float-depth))))

;;; Could change to ',' or something special.
(defun delimiterp (c)
  (or (char= c #\Space)
      (char= c #\Tab)))

(defun split-per-space (str &key (delimiterp #'delimiterp))
  "Split string with specified delimiter."
  (loop for point = (position-if-not delimiterp str)
     then (position-if-not delimiterp str :start (1+ end))
     for end = (and point (position-if delimiterp str :start point))
     when point collect (subseq str point end)
     while end))

(defun make-array-per-space (str)
  "Make an array from a string with words seperated by delimiters."
  (let ((nums nil))
    (dolist (num (split-per-space str))
      (push (parse-decimal num) nums))
    (list-array (reverse nums))))

(defun make-matrix-per-space (str-lst)
  "Make a matrix from a list of strings with words seperated by delimiters."
  (let ((nums-m nil)
	(nums nil))
    (dolist (str str-lst)
      (progn
	(setf nums nil)
	(dolist (num (split-per-space str))
	  (push (parse-decimal num) nums))
	(progn
	  (setf nums (reverse nums))
	  (push nums nums-m))))
    (list-matrix (reverse nums-m))))

;;; Will be useful to implement the Baum-Welch algorithm.
(defun get-gamma (hstate tms alpha beta)
  (declare (type fixnum tms)
	   (type matrix alpha)
	   (type matrix beta))
  (let ((denom 0.0)
	(gamma (make-array (list tms (get-states hstate)) :element-type 'float)))
    (do ((tm 0 (1+ tm)))
	((>= tm (- tms 1)))
      (progn
	(setf denom 0.0)
	(do ((i 0 (1+ i)))
	    ((>= i (get-states hstate)))
	  (progn
	    (setf (aref gamma tm i) (* (aref alpha tm i) (aref beta tm i)))
	    (incf denom (aref gamma tm i))))
	(do ((j 0 (1+ j)))
	    ((>= j (get-states hstate)))
	  (setf (aref gamma tm j) (/ (aref gamma tm j) denom)))))))

(defun get-xi (hstate tms o alpha beta)
  (declare (type fixnum tms)
	   (type array o)
	   (type matrix alpha)
	   (type matrix beta))
  (let ((sum 0.0)
	(xi (make-array (list tms (get-states hstate) (get-states hstate)) :element-type 'float)))
    (do ((tm 0 (1+ tm)))
	((>= tm (- tms 1)))
      (progn
	(setf sum 0.0)
	(do ((i 0 (1+ i)))
	    ((>= i (get-states hstate)))
	  (do ((j 0 (1+ j)))
	      ((>= j (get-states hstate)))
	    (progn
	      (setf (aref xi tm i j) (* (aref alpha tm i)
					(aref beta (+ tm 1) j)
					(aref (get-transition hstate) i j)
					(aref (get-measure hstate) j (aref o (+ tm 1)))))
	      (incf sum (aref xi tm i j)))))
	(do ((k 0 (1+ k)))
	    ((>= k (get-states hstate)))
	  (do ((l 0 (1+ l)))
	      ((>= l (get-states hstate)))
	    (setf (aref xi tm k l) (/ (aref xi tm k l) sum))))))))

;;; APIs
(defun init-hmm (s o i tr m)
  "Initiate an hmm-state with given arguments."
  (declare (type fixnum s)
	   (type fixnum o)
	   (type vector i)
	   (type matrix tr)
	   (type matrix m))
  (make-hmm-state :states s :measures o :initial i :transition-matrix tr
		  :measure-matrix m))

(defun init-from-file (file)
  "Initiate a brand new hmm-state from a input file."
  (with-open-file (in file :direction :input)
    (let ((hstate (make-hmm-state)))
      (do ((line (read-line in nil)
		 (read-line in nil)))
	  ((null line))
	(cond ((string= "S" line :end2 1)
	       ;; states
	       (setf (get-states hstate) (parse-integer line :start 4)))
	      ((string= "O" line :end2 1)
	       ;; observations
	       (setf (get-observations hstate) (parse-integer line :start 4)))
	      ((string= "I" line :end2 1)
	       ;; initial distributions
	       (setf line (read-line in nil)
		     (get-initial hstate) (make-array-per-space line)))
	      ((string= "T" line :end2 1)
	       ;; transition matrix
	       (let ((str-lst nil))
		 (progn
		   (do ((str-line (read-line in nil)
				  (read-line in nil)))
		       ((string= #\Space str-line :end2 1))
		     (push str-line str-lst))
		   (setf (get-transition hstate) (make-matrix-per-space (reverse str-lst))))))
	      ((string= "M" line :end2 1)
	       ;; measure matrix
	       (let ((str-lst nil))
		 (progn
		   (do ((str-line (read-line in nil)
				  (read-line in nil)))
		       ((null str-line))
		     (push str-line str-lst))
		   (setf (get-measure hstate) (make-matrix-per-space (reverse str-lst))))))
	      ((string= #\Space line))
	      ;; jump over space
	      (t
	       (error "broken HMM file."))))
    hstate)))

(defun print-to-file (hstate file)
  "Store an hmm-state information into a file."
  (with-open-file (out file
		       :direction :output
		       :if-exists :supersede)
    (progn
      (format out "S = ~d~%" (get-states hstate))
      (format out " ~%")
      (format out "O = ~d~%" (get-observations hstate))
      (format out " ~%")
      (format out "I =~%")
      (pprint-array out (get-initial hstate))
      (format out "~%")
      (format out " ~%")
      (format out "T =~%")
      (pprint-matrix out (get-transition hstate))
      (format out " ~%")
      (format out "M =~%")
      (pprint-matrix out (get-measure hstate)))))

(defun forward (hstate tms o)
  "Forward algorithm."
  (declare (type array o)
	   (type fixnum tms))
  (let ((alpha (make-array (list tms (get-states hstate)) :element-type 'float))
	(sum 0.0)
	(prob 0.0))
    (progn
      (do ((i 0 (1+ i)))
	  ((>= i (get-states hstate)))
	(setf (aref alpha 0 i) (* (aref (get-initial hstate) i)
				  (aref (get-measure hstate) i (aref o 0)))))
      (do ((tm 0 (1+ tm)))
	  ((>= tm (- tms 1)))
	(do ((j 0 (1+ j)))
	    ((>= j (get-states hstate)))
	  (progn
	    (setf sum 0.0)
	    (do ((k 0 (1+ k)))
		((>= k (get-states hstate)))
	      (incf sum (* (aref alpha tm k)
			   (aref (get-transition hstate) k j))))
	    (setf (aref alpha (+ 1 tm) j) (* sum (aref (get-measure hstate) j (aref o (+ 1 tm))))))))
      (do ((l 0 (1+ l)))
	  ((>= l (get-states hstate)))
	(incf prob (aref alpha (- tms 1) l))))
    (values alpha prob)))

(defun backward (hstate tms o)
  "Backward algorithm."
  (declare (type array o)
	   (type fixnum tms))
  (let ((beta (make-array (list tms (get-states hstate)) :element-type 'float))
	(sum 0.0)
	(prob 0.0))
    (progn
      (do ((i 0 (1+ i)))
	  ((>= i (get-states hstate)))
	(setf (aref beta (- tms 1) i) 1.0))
      (do ((tm (- tms 2) (1- tm)))
	  ((< tm 0))
	(do ((k 0 (1+ k)))
	    ((>= k (get-states hstate)))
	  (progn
	    (setf sum 0.0)
	    (do ((j 0 (1+ j)))
		((>= j (get-states hstate)))
	      (incf sum (* (aref (get-transition hstate) k j)
			   (aref (get-measure hstate) j (aref o (+ 1 tm)))
			   (aref beta (+ 1 tm) j))))
	    (setf (aref beta tm k) sum))))
      (do ((l 0 (1+ l)))
	  ((>= l (get-states hstate)))
	(incf prob (aref beta 0 l))))
    (values beta prob)))

(defun viterbi (hstate tms o)
  "Viterbi algorithm."
  (declare (type array o)
	   (type fixnum tms))
  (let ((delta (make-array (list tms (get-states hstate)) :element-type 'float))
	(psi (make-array (list tms (get-states hstate)) :element-type 'integer))
	(q (make-array tms :element-type 'integer))
	(max-value 0.0)
	(max-value-indice 0)
	(value 0.0)
	(prob 0.0))
    (progn
      (do ((i 0 (1+ i)))
	  ((>= i (get-states hstate)))
	(setf (aref delta 0 i) (* (aref (get-initial hstate) i)
				  (aref (get-measure hstate) i (aref o 0)))
	      (aref psi 0 i) 0))
      (do ((tm 1 (1+ tm)))
	  ((>= tm tms))
        (do ((j 0 (1+ j)))
	    ((>= j (get-states hstate)))
	  (progn
	    (setf max-value 0.0
		  max-value-indice 1)
	    (do ((k 0 (1+ k)))
		((>= k (get-states hstate)))
	      (progn
		(setf value (* (aref delta (- tm 1) k)
			       (aref (get-transition hstate) k j)))
		(if (> value max-value)
		    (setf max-value value
			  max-value-indice k))))
	    (setf (aref delta tm j) (* max-value (aref (get-measure hstate) j (aref o tm)))
		  (aref psi tm j) max-value-indice))))
      (setf (aref q (- tms 1)) 1)
      (do ((l 0 (1+ l)))
	  ((>= l (get-states hstate)))
	(if (> (aref delta (- tms 1) l) prob)
	    (setf prob (aref delta (- tms 1) l)
		  (aref q (- tms 1)) l)))
      (do ((tm (- tms 2) (1- tm)))
	  ((< tm 0))
	(setf (aref q tm) (aref psi (+ tm 1) (aref q (+ tm 1))))))
    (values q prob)))
