;; Produces output for gsl_const header files using GNU Calc.
;;
;; Generate output with
;;
;;   emacs -batch -l const.el -f run 
;;

(setq calc-display-working-message t) ;; display short working messages
(setq calc-float-format '(sci 20))
(calc-eval "")
(load-library "calc/calc-units.el")
(calc-extensions)

(setq  gsl-constants
  '(("c" "SPEED_OF_LIGHT")
    ("au" "ASTRONOMICAL_UNIT")
    ("Grav" "GRAVITATIONAL_CONSTANT"))
  )

;;; work around bug in calc 2.02f
(defun math-extract-units (expr)
  (if (memq (car-safe expr) '(* /))
      (cons (car expr)
	    (mapcar 'math-extract-units (cdr expr)))
    (if (math-units-in-expr-p expr nil) expr 1))
)



;;  (nth 1 (assq 'cgs math-standard-units-systems))




(defun fn (prefix system expr name)
  (let* ((x (calc-eval expr 'raw))
         (y (math-to-standard-units x system))
         (z (math-simplify-units y))
         (quantity (calc-eval (math-remove-units z)))
         (units (calc-eval (math-extract-units z)))
         )
    (print x)
    (print y)
    (print z)
    (print (math-extract-units z))
    (print quantity)
    (print units)
    (princ (format "#define %s_%s (%s) /* %s */\n" prefix name quantity units))
    )
  )

(setq cgs (nth 1 (assq 'cgs math-standard-units-systems)))


(defun run ()
  (mapcar (lambda (x) (apply 'fn "GSL_CONST" cgs x)) gsl-constants)
)

