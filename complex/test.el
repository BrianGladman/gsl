;; Produces test output for complex functions
;; using GNU Calc.
;;
;; Generate output with
;;
;;   emacs -batch -l test.el -f test-all > results.h
;;
;; Note: this takes a long time to run

;; set Calc to use radians mode, turn off symbolic
;; evaluation and use a reasonable precision.

(setq max-lisp-eval-depth 4000)
(setq max-specpdl-size 2400)
(setq calc-internal-prec 64) ;; probably unnecessarily high, but found a few discrepancies at prec=20
(setq calc-angle-mode 'rad)
(setq calc-float-format '(sci 20))
;;(setq calc-full-float-format '(sci 0))

(setq z '("0.0" "1.0e-8" "0.1" "0.9" "1.0" "10.0" "-1.0e-8"
          "-0.1" "-0.9" "-1.0" "-10.0"))

;; faster
(setq z '("0.0" "1.0e-8" "0.9" "10.0" "-1.0e-8" "-0.9" "-10.0"))

(setq f '(;;"arg" "abs" 
          ;;"sqrt" "log" "log10" "exp" 
          ;; "sin" "cos" "tan" 
          "arcsin" "arccos" "arctan" 
          "sinh" "cosh" "tanh" "arcsinh" "arccosh" "arctanh"
          ))

(defun test (function arg1 arg2)
  (let* ((v (concat function "((" arg1 "," arg2 "))"))
        (result (calc-eval v)))
    (if (string-match "(\\(.*\\),\\(.*\\))" result)
        (setq result (replace-match "\\1,\\2" nil nil result)))
    (if (string-match "^\\([^,]*\\)$" result)
        (setq result (replace-match "\\1, 0.0" nil nil result)))
    (if (not (string-match "(" result))  ;; skip any unsimplified results
        (princ (format "  {FN (%s), %s, %s, %s},\n" function arg1 arg2 result))
      )))

;;(test "sin" "10" "0")

;; loop over all possible combinations of a,b,c

(defun loop (a b c)
  (while a
    (progn 
      (let ((b1 b))
        (while b1
          (progn
            (let ((c1 c))
              (while c1
                (progn 
                  (test (car a) (car b1) (car c1))
                  (setq c1 (cdr c1))
                  )
                )
              )
            (setq b1 (cdr b1))
            )
          )
        )
      (setq a (cdr a))
      )
    )
  )

(defun test-all ()
  (loop f z z)
)

;;(test-all)



