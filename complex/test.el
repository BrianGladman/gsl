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

(setq z '(0.0 
          1.0e-8 
          0.1 
          0.5 
          0.9 
          1.0
          2.0
          10.0 
          3.14159265358979323846264338328   ;; pi
          1.57079632679489661923132169164   ;; pi/2
          0.78539816339744830966156608458   ;; pi/4
          4.7123889803846898576939650749193   ;; 3pi/2
          1.77245385090551602729816748334   ;; sqrt(pi)
          1.14472988584940017414342735135   ;; ln(pi) 
          ))

(defun reflections (a b)
  (let ((a (float a)) (b (float b)))
    (list
     (format "(%.16g,%.16g)" a b)
     (format "(%.16g,%.16g)" a (- b))
     (format "(%.16g,%.16g)" (- a) b)
     (format "(%.16g,%.16g)" (- a) (- b))
     (format "(%.16g,%.16g)" b a)
     (format "(%.16g,%.16g)" b (- a))
     (format "(%.16g,%.16g)" (- b) a)
     (format "(%.16g,%.16g)" (- b) (- a)))))

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

;;(trim (permute 'reflections (list 0 1) (list 0 1)))

;; faster
;;(setq z '("0.0" "1.0e-8" "0.9" "10.0" "-1.0e-8" "-0.9" "-10.0"))

(setq f '(;;"arg" "abs" 
          "sqrt" "log" "log10" "exp" 
          "sin" "cos" "tan" 
          "arcsin" "arccos" "arctan" 
          "sinh" "cosh" "tanh" "arcsinh" "arccosh" "arctanh"
          ))

(defun test (function arg)
  (let* ((v (concat "clean(" function "(" arg "),60)"))
        (result (calc-eval v)))
    (if (string-match "(\\(.*\\),\\(.*\\))" result)
        (setq result (replace-match "\\1,\\2" nil nil result)))
    (if (string-match "^\\([^,]*\\)$" result)
        (setq result (replace-match "\\1, 0.0" nil nil result)))
    (if (not (string-match "(" result))  ;; skip any unsimplified results
        (princ (format "  {FN (%s), ARG%s, RES(%s)},\n" function arg result))
      )))

;;(test "sin" "10" "0")

;; loop over all possible combinations of a,b,c

(defun loop (a b)
  (while a
    (progn 
      (let ((b1 b))
        (while b1
          (progn
            (message "%s(%s)" (car a) (car b1))
            (test (car a) (car b1))
            (setq b1 (cdr b1))
            )
          )
        )
      (setq a (cdr a))
      )
    )
  )

(defun test-all ()
  (loop f (trim (permute 'reflections z z)))
)

;;(test-all)




