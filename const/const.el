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
       '(("fsc"           "FINE_STRUCTURE")
         ("Nav"           "AVAGADRO")
         )
       )

(setq  gsl-constants
       '(("c"             "SPEED_OF_LIGHT")
         ("Grav"          "GRAVITATIONAL_CONSTANT")
         ("h"             "PLANCKS_CONSTANT_H")
         ("hbar"          "PLANCKS_CONSTANT_HBAR")
         ("mu0"           "VACUUM_PERMEABILITY")

         ("au"            "ASTRONOMICAL_UNIT")
         ("float(lyr)"    "LIGHT_YEAR")
         ("pc"            "PARSEC")

         ("ga"            "GRAV_ACCEL")

         ("ev"            "ELECTRON_VOLT")
         ("me"            "MASS_ELECTRON")
         ("mu"            "MASS_MUON")
         ("mp"            "MASS_PROTON")
         ("mn"            "MASS_NEUTRON")

         ("Ryd"           "RYDBERG")
         ("k"             "BOLTZMANN")
         ("muB"           "BOHR_MAGNETON")
         ("muN"           "NUCLEAR_MAGNETON")
         ("mue"           "ELECTRON_MAGNETIC_MOMENT")
         ("mup"           "PROTON_MAGNETIC_MOMENT")

         ("R0"            "MOLAR_GAS")
         ("V0"            "STANDARD_GAS_VOLUME")

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
         ("gal"           "US_GALLON")
         ("qt"            "QUART")
         ("pt"            "PINT")
         ("cup"           "CUP")
         ("ozfl"          "FLUID_OUNCE")
         ("tbsp"          "TABLESPOON")
         ("tsp"           "TEASPOON")
         ("galC"          "CANADIAN_GALLON")
         ("galUK"         "UK_GALLON")
         
         ("mph"           "MILES_PER_HOUR")
         ("kph"           "KILOMETERS_PER_HOUR")
         ("knot"          "KNOT")

         ("lb"            "POUND_MASS")
         ("oz"            "OUNCE_MASS")
         ("ton"           "TON")
         ("t"             "METRIC_TON")
         ("tonUK"         "UK_TON")
         ("ozt"           "TROY_OUNCE")
         ("ct"            "CARAT")
         ("amu"           "UNIFIED_ATOMIC_MASS")

         ("gf"            "GRAM_FORCE")
         ("lbf"           "POUND_FORCE")
         ("kip"           "KILOPOUND_FORCE")
         ("pdl"           "POUNDAL")

         ("cal"           "CALORIE")
         ("btu"           "BTU")
         ("ev"            "ELECTRON_VOLT")
         ("therm"         "THERM")

         ("hp"            "HORSEPOWER")
         
         ("bar"           "BAR")
         ("atm"           "STD_ATMOSPHERE")
         ("torr"          "TORR")
         ("mHg"           "METER_OF_MERCURY")
         ("inHg"          "INCH_OF_MERCURY")
         ("inH2O"         "INCH_OF_WATER")
         ("psi"           "PSI")

         ("P"             "POISE")
         ("St"            "STOKES")
         
         ("Fdy"           "FARADAY")
         ("e"             "ELECTRON_CHARGE")
         ("G"             "GAUSS")

         ("sb"            "STILB")
         ("lm"            "LUMEN")
         ("lx"            "LUX")
         ("ph"            "PHOT")
         ("fc"            "FOOTCANDLE")
         ("lam"           "LAMBERT")
         ("flam"          "FOOTLAMBERT")
         
         ("Ci"            "CURIE")
         ("R"             "ROENTGEN")
         ("rd"            "RAD")


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

