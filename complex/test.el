;; Produces test output for complex functions using GNU Calc.
;;
;; Generate output with
;;
;;   emacs -batch -l test.el -f test-all > results.h
;;
;; Note: this takes a long time to run

;; set Calc to use radians mode, turn off symbolic evaluation and use
;; a reasonable precision.

(setq calc-display-working-message t) ;; display short working messages
(setq max-lisp-eval-depth 4000)
(setq max-specpdl-size 2400)
(setq calc-internal-prec 64) ;; probably unnecessarily high, but found
			     ;; a few discrepancies at prec=20
(setq calc-angle-mode 'rad)
(setq calc-float-format '(sci 20))
;;(setq calc-full-float-format '(sci 0))

(defun reflections (a b)
  (let ((a (float a)) (b (float b)))
    (list
     (list a b)
     (list a (- b))
     (list (- a) b)
     (list (- a) (- b))
     (list b a)
     (list b (- a))
     (list (- b) a)
     (list (- b) (- a)))))

(defun permute (fn a b)
  (let ((a1 a) (result nil))
    (while a1
      (progn 
        (let ((b1 b))
          (while b1
            (progn 
              (setq result (append result (funcall fn (car a1) (car b1))))
              (setq b1 (cdr b1))
              )
            )
          )
        (setq a1 (cdr a1))
        )
      )
    result
    )
  )

(defun trim (a)
  (let ((result nil))
    (while a
      (setq result (cons (car a) result))
      (setq a (delete (car a) a))
      )
    (reverse result))
)

(defun combine (a b)
  (trim (permute 'reflections a b)))

(defun evaltest (function arg)
  (let* ((x (nth 0 arg))
         (y (nth 1 arg))
         (z (format "(%.60e,%.60e)" x y))
         (v (concat "clean(" function "(" z "),60)"))
        (result (calc-eval v)))
    ;;(message "%s (%g %g)" function x y)
    (if (string-match "(\\(.*\\),\\(.*\\))" result)
        (setq result (replace-match "\\1,\\2" nil nil result)))
    (if (string-match "^\\([^,]*\\)$" result)
        (setq result (replace-match "\\1, 0.0" nil nil result)))
    (if (not (string-match "(" result))  ;; skip any unsimplified results
        (princ (format "  {FN (%s), ARG(%.18e,%.18e), RES(%s)},\n" function x y result))
      )))

;;(evaltest "sin" "10" "0")

;; loop over all possible combinations of a,b,c

(defun test (a b)
  (let ((b1 b))
    (while b1
      (progn
        (evaltest a (car b1))
        (setq b1 (cdr b1))
        )
      )
    )
  )
  
;;

(setq pi 3.14159265358979323846264338328);
(setq flteps 1.1920928955078125e-07);
(setq eps (sqrt flteps));

(setq edge (list flteps))
(setq simple (list 0.0 0.5  1.0  2.0  10.0))
(setq real (list 0.0 0.125 0.5 0.75  1.0  2.0  10.0))
(setq trig (list (+ pi eps)
                 (+ (* 0.5 pi) eps)
                 (+ (* 0.25 pi) eps)
                 (+ (* 0.75 pi) eps)
                 (+ (sqrt pi) eps)
                 (+ (log pi) eps)
                 (- pi eps)
                 (- (* 0.5 pi) eps)
                 (- (* 0.25 pi) eps)
                 (- (* 0.75 pi) eps)
                 (- (sqrt pi) eps)
                 (- (log pi) eps)))

(setq reals (append (combine simple real) 
                    (combine edge (append edge simple) )))

(setq trigs (append reals 
                    (combine simple trig)
                    (combine edge trig)))

(defun test-all ()

;  (test "sqrt" reals)

;  (test "log" reals)
;  (test "log10" reals)
;  (test "exp" trigs)

;  (test "sin" trigs)
;  (test "cos" trigs)
;  (test "tan" trigs)

;  (test "arcsin" trigs)
;  (test "arccos" trigs)
;  (test "arctan" trigs)

;  (test "sinh" trigs)
;  (test "cosh" trigs)
;  (test "tanh" trigs)

;  (test "arcsinh" trigs)
;  (test "arccosh" trigs)
  (test "arctanh" trigs)

  (test "cosec" trigs)
  (test "sec" trigs)
  (test "cot" trigs)

  (test "cosech" trigs)
  (test "sech" trigs)
  (test "coth" trigs)

  (test "arccosec" trigs)
  (test "arcsec" trigs)
  (test "arccot" trigs)

  (test "arccosech" trigs)
  (test "arcsech" trigs)
  (test "arccoth" trigs)

)

;;(test-all)




