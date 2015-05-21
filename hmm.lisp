;;;; Common Lisp Hidden Markov Model Algorithms Package

(defpackage :hmm-algorithms
  (:nicknames :hmm :cl-hmm)
  (:use :cl)
  (:export :init-from-file
	   :print-to-file
	   :forward
	   :backward
	   :viterbi))

(in-package :hmm-algorithms)

;;; Matrix type, in fact a 2-dimensional array.
(deftype matrix (&optional type x y)
  "Transition square matrix and measure rectangle matrix."
  `(array ,type (,x ,y)))

;;; In order to pretty print an array, implement it.
(defun pprint-array (stream array)
  "Pretty print routine to print a nested vector."
  (loop
     :with first-time = t
     :for x :across array
     :unless first-time
     :do (write-char #\Space stream) :end
     :do (princ x stream)
     (setf first-time nil)))

;;; And so with a matrix...
(defun pprint-matrix (stream matrix)
  (loop for i below (car (array-dimensions matrix)) do
       (loop for j below (cadr (array-dimensions matrix)) do
	    (let ((cell (aref matrix i j)))
	      (format stream "~a " cell)))
       (format stream "%")))

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
(defmacro get-states (hstate)
  `(hmm-state-states ,hstate))

(defmacro get-observations (hstate)
  `(hmm-state-measures ,hstate))

(defmacro get-initial-distributions (hstate)
  `(hmm-state-initial ,hstate))

(defmacro get-transition-matrix (hstate)
  `(hmm-state-transition-matrix ,hstate))

(defmacro get-measure-matrix (hstate)
  `(hmm-state-measure-matrix ,hstate))

;;; Error handling routines.
(define-condition input-format-error (error)
  ((text :initarg :text :reader text)))

;;; Helpers, to generate a new HMM state from file.
(defun list-array (lst)
  (make-array (length lst) :initial-contents lst))

(defun list-matrix (lst)
  (make-array (list (length lst)
		    (length (first lst)))
	      :initial-contents lst))

(defun delimiterp (c)
  (or (char= c #\Space)
      (char= c #\Tab)))

(defun split-per-space (str &key (delimiterp #'delimiterp))
  (loop :for point = (position-if-not delimiterp str)
     :then (position-if-not delimiterp str :start (1+ end))
     :for end = (and point (position-if delimiterp str :start point))
     :when point :collect (subseq str point end)
     :while end))

(defun make-array-per-space (str)
  (let ((nums nil))
    (dolist (num (split-per-space str))
      (push (parse-integer num) nums))
    (list-array (reverse nums))))

(defun make-matrix-per-space (str-lst)
  (let ((nums-m nil)
	(nums nil))
    (dolist (str str-lst)
      (progn
	(setf nums nil)
	(dolist (num (split-per-space str))
	  (push (parse-integer num) nums))
	(progn
	  (setf nums (reverse nums))
	  (push nums nums-m))))
    (list-matrix (reverse nums-m))))
  
(defun get-gamma (hstate times alpha beta)
  (declare (type fixnum times)
	   (type matrix alpha)
	   (type matrix beta))
  (let ((denom 0.0d0)
	(gamma (make-array '(times (get-states hstate)) :element-type 'double-float)))
    (do ((tm 0 (1+ tm)))
	((>= tm (- times 1)))
      (progn
	(setf denom 0.0d0)
	(do ((i 0 (1+ i)))
	    ((>= i (get-states hstate)))
	  (progn
	    (setf (aref gamma tm i) (* (aref alpha tm i) (aref beta tm i)))
	    (incf denom (aref gamma tm i))))
	(do ((j 0 (1+ j)))
	    ((>= j (get-states hstate)))
	  (setf (aref gamma tm j) (/ (aref gamma tm j) denom)))))))

(defun get-xi (hstate times o alpha beta)
  (declare (type fixnum times)
	   (type array o)
	   (type matrix alpha)
	   (type matrix beta))
  (let ((sum 0.0d0)
	(xi (make-array '(times (get-states hstate) (get-states hstate)) :element-type 'double-float)))
    (do ((tm 0 (1+ tm)))
	((>= tm (- times 1)))
      (progn
	(setf sum 0.0d0)
	(do ((i 0 (1+ i)))
	    ((>= i (get-states hstate)))
	  (do ((j 0 (1+ j)))
	      ((>= j (get-states hstate)))
	    (progn
	      (setf (aref xi tm i j) (* (aref alpha tm i)
					(aref beta (+ tm 1) j)
					(aref (get-transition-matrix hstate) i j)
					(aref (get-measure-matrix hstate) j (aref o (+ tm 1)))))
	      (incf sum (aref xi tm i j)))))
	(do ((k 0 (1+ k)))
	    ((>= k (get-states hstate)))
	  (do ((l 0 (1+ l)))
	      ((>= l (get-states hstate)))
	    (setf (aref xi tm k l) (/ (aref xi tm k l) sum))))))))

;;; APIs
(defun init-from-file (file)
  "Initiate a brand new hmm state from a input file."
  (with-open-file (in file :direction :input)
    (let ((hstate (make-hmm-state)))
      (do ((line (read-line in nil)
		 (read-line in nil)))
	  ((null line))
	(progn
	  (setf line (string-trim '(#\Space #\Tab) line))
	  (cond ((string= "S" line :end2 1)
		 (setf (get-states hstate) (parse-integer line :start 4)))
		((string= "O" line :end2 1)
		 (setf (get-observations hstate) (parse-integer line :start 4)))
		((string= "I" line :end2 1)
		 (progn
		   (setf line (string-trim '(#\Space #\Tab) (read-line in nil)))
		   (if (null line)
		       (error 'input-format-error :text "broken HMM file.")
		       (setf (get-initial-distributions hstate)
			     (make-array-per-space line)))))
		((string= "T" line :end2 1)
		 (let ((str-lst nil))
		   (progn
		     (do ((str-line (read-line in nil)
				    (read-line in nil)))
			 ((or (string= ">" str-line :end2 1)
			      (null str-line)))
		       (push line str-lst))
		     (setf (get-transition-matrix hstate) (make-matrix-per-space (reverse str-lst))))))
		((string= "M" line :end2 1)
		 (let ((str-lst nil))
		   (progn
		     (do ((str-line (read-line in nil)
				    (read-line in nil)))
			 ((or (string= ">" str-line :end2 1)
			      (null str-line)))
		       (push line str-lst))
		     (setf (get-measure-matrix hstate) (make-matrix-per-space (reverse str-lst))))))	  		   
		(t
		 (error 'input-format-error
			:text "broken HMM file.")))))
      hstate)))

(defun print-to-file (hstate file)
  "Store a hmm state information into a file."
  (with-open-file (out file
		       :direction :output
		       :if-exists :supersede)
    (progn
      (format out "S = ~d~%" (get-states hstate))
      (format out "O = ~d~%" (get-observations hstate))
      (format out "I = <~%")
      (pprint-array out (get-initial-distributions hstate))
      (format out ">~%")
      (format out "T = <~%")
      (pprint-matrix out (get-transition-matrix hstate))
      (format out ">~%")
      (format out "M = <~%")
      (pprint-matrix out (get-measure-matrix hstate))
      (format out ">~%"))))

(defun forward (hstate times o)
  "Forward algorithm."
  (declare (type fixnum times)
	   (type array o))
  (let ((alpha (make-array '(times (get-states hstate)) :element-type 'double-float))
	(sum 0.0d0)
	(prob 0.0d0))
    (progn
      (do ((i 0 (1+ i)))
	  ((>= i (get-states hstate)))
	(setf (aref alpha 0 i) (* (aref (get-initial-distributions hstate) i)
				  (aref (get-measure-matrix hstate) i (aref o 0)))))
      (do ((tm 0 (1+ tm)))
	  ((>= tm (- times 1)))
	(do ((j 0 (1+ j)))
	    ((>= j (get-states hstate)))
	  (progn
	    (do ((k 0 (1+ k)))
		((>= k (get-states hstate)))
	      (incf sum (* (aref alpha tm k)
			   (aref (get-transition-matrix hstate) k j))))
	    (setf (aref alpha tm j) (* sum (aref (get-measure-matrix hstate) j (aref o tm)))))))
      (do ((l 0 (1+ l)))
	  ((>= l (get-states hstate)))
	(incf prob (aref alpha (- times 1) l))))
    (values alpha prob)))

(defun backward (hstate times o)
  "Backward algorithm."
  (declare (type fixnum times)
	   (type array o))
  (let ((beta (make-array '(times (get-states hstate)) :element-type 'double-float))
	(sum 0.0d0)
	(prob 0.0d0))
    (progn
      (do ((i 0 (1+ i)))
	  ((>= i (get-states hstate)))
	(setf (aref beta (- times 1) i) 1.0))
      (do ((tm (- times 2) (1- tm)))
	  ((< tm 0))
	(do ((k 0 (1+ k)))
	    ((>= k (get-states hstate)))
	  (progn
	    (do ((j 0 (1+ j)))
		((>= j (get-states hstate)))
	      (incf sum (* (aref (get-transition-matrix hstate) k j)
			   (aref (get-measure-matrix hstate) j (aref o tm))
			   (aref beta tm j))))
	    (setf (aref beta (- tm 1) k) sum))))
      (do ((l 0 (1+ l)))
	  ((>= l (get-states hstate)))
	(incf prob (aref beta 0 l))))
    (values beta prob)))

(defun viterbi (hstate times o)
  "Viterbi algorithm."
  (declare (type fixnum times)
	   (type array o))
  (let ((delta (make-array '(times (get-states hstate)) :element-type 'double-float))
	(psi (make-array '(times (get-states hstate)) :element-type 'integer))
	(q (make-array times :element-type 'integer))
	(max-value 0.0d0)
	(max-value-indice 0)
	(value 0.0d0)
	(prob 0.0d0))
    (progn
      (do ((i 0 (1+ i)))
	  ((>= i (get-states hstate)))
	(progn
	  (setf (aref delta 0 i) (* (aref (get-initial-distributions hstate) i)
				    (aref (get-measure-matrix hstate) i (aref o 0))))
	  (setf (aref psi 0 i) 0)))
      (do ((tm 1 (1+ tm)))
	  ((>= tm times))
        (do ((j 0 (1+ j)))
	    ((>= j (get-states hstate)))
	  (progn
	    (setf max-value 0.0d0
		  max-value-indice 1)
	    (do ((k 0 (1+ k)))
		((>= k (get-states hstate)))
	      (progn
		(setf value (* (aref delta (- tm 1) k)
			       (aref (get-transition-matrix hstate) k j)))
		(if (> value max-value)
		    (setf max-value value
			  max-value-indice k))))
	    (setf (aref delta tm j) (* max-value (aref (get-measure-matrix hstate) j (aref o tm)))
		  (aref psi tm j) max-value-indice))))
      (setf (aref q (- times 1)) 1)
      (do ((l 0 (1+ l)))
	  ((>= l (get-states hstate)))
	(if (> (aref delta (- times 1) l) prob)
	    (setf prob (aref delta (- times 1) l)
		  (aref q (- times 1)) l)))
      (do ((tm (- times 2) (1- tm)))
	  ((< tm 0))
	(setf (aref q tm) (aref psi (+ tm 1) (aref q (+ tm 1))))))
    (values q prob)))