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

(setq  gsl-dimensionless-constants
       '(("fsc"           "FINE_STRUCTURE_ALPHA")
         ("Nav"           "AVAGADRO")
         )
       )

(setq  gsl-constants
       '(("c"             "SPEED_OF_LIGHT")
         ("Grav"          "GRAVITATIONAL_CONSTANT")
         ("h"             "PLANCKS_CONSTANT_H")
         ("hbar"          "PLANCKS_CONSTANT_HBAR")

         ("au"            "ASTRONOMICAL_UNIT")
         ("float(lyr)"    "LIGHT_YEAR")
         ("pc"            "PARSEC")

         ("ga"            "GRAV_ACCEL")

         ("ev"            "ELECTRON_VOLT")
         ("me"            "MASS_ELECTRON")
         ("mu"            "MASS_MUON")
         ("mp"            "MASS_PROTON")
         ("mn"            "MASS_NEUTRON")

         ("min"           "MINUTE")
         ("hr"            "HOUR")
         ("day"           "DAY")
         ("wk"            "WEEK")

         ("in"            "INCH")
         ("ft"            "FOOT")
         ("yd"            "YARD")
         ("mi"            "MILE")
         ("nmi"           "NAUTICAL_MILE")
         ("fath"          "FATHOM")

         ("mil"           "MIL")
         ("pt"            "POINT")
         ("tpt"           "TEXPOINT")

         ("u"             "MICRON")
         ("Ang"           "ANGSTROM")

         ("hect"          "HECTARE")
         ("acre"          "ACRE")
         ("barn"          "BARN")

         ("l"             "LITER")
         ("gal"           "USGALLON")
         ("qt"            "QUART")
         ("pt"            "PINT")
         ("cup"           "CUP")
         ("ozfl"          "FLUID_OUNCE")
         ("tbsp"          "TABLESPOON")
         ("tsp"           "TEASPOON")
         ("galC"          "CANADIANGALLON")
         ("galUK"         "UKGALLON")
         
         ("mph"           "MILES_PER_HOUR")
         ("kph"           "KILOMETERS_PER_HOUR")
         ("knot"          "KNOT")

         ("amu"            "UNIFIED_ATOMIC_MASS")

         )
       )

;;; work around bug in calc 2.02f
(defun math-extract-units (expr)
  (if (memq (car-safe expr) '(* /))
      (cons (car expr)
	    (mapcar 'math-extract-units (cdr expr)))
    (if (math-units-in-expr-p expr nil) expr 1))
)

(defun fn (prefix system expr name)
  (let* ((x (calc-eval expr 'raw))
         (y (math-to-standard-units x system))
         (z (math-simplify-units y))
         (quantity (calc-eval (math-remove-units z)))
         (units (calc-eval (math-extract-units z)))
         )
    ;;(print x)
    ;;(print y)
    ;;(print z)
    ;;(print (math-extract-units z))
    ;;(print quantity)
    ;;(print units)
    (princ (format "#define %s_%s (%s) /* %s */\n" prefix name quantity units))
    )
  )

(setq cgs (nth 1 (assq 'cgs math-standard-units-systems)))
(setq mks (nth 1 (assq 'mks math-standard-units-systems)))

(defun display (prefix system constants)
  (princ (format "#ifdef __%s__\n" prefix))
  (princ (format "#define __%s__\n\n" prefix))
  (mapcar (lambda (x) (apply 'fn prefix system x)) constants)
  (princ (format "\n#endif /* __%s__ */\n" prefix))
)

(defun run-cgs ()
  (display "GSL_CONST_CGS" cgs gsl-constants)
)

(defun run-mks ()
  (display "GSL_CONST_MKS" mks gsl-constants)
)

(defun run-num ()
  (display "GSL_CONST_NUM" mks gsl-dimensionless-constants)
)
