;; Produces output for gsl_const header files using GNU Calc.
;;
;; Generate output with
;;
;;   emacs -batch -l test.el -f run 
;;

(setq calc-display-working-message t) ;; display short working messages
(calc-eval "")
(calc-extensions)
(setq calc-float-format '(sci 20))

(setq  gsl-constants
  '(("c" "SPEED_OF_LIGHT")
    ("au" "ASTRONOMICAL_UNIT"))
  )

(defun fn (expr name)
  (let ((value (calc-eval 
                (math-remove-units 
                 (math-to-standard-units (calc-eval expr 'raw) nil)))))
    (princ (format "#define GSL_CONST_%s (%s)\n" name value))
    )
  )

(defun run ()
  (mapcar (lambda (x) (apply 'fn x)) gsl-constants)
)

;;(apply 'fn (nth 1 gsl-constants))

;;(fn "au" "x")