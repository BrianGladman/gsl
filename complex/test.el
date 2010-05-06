;; Produces test output for complex functions using GNU Calc.
;;
;; Generate output with
;;
;;   emacs -batch -l test.el -f test1 | sed  's/00\+e/0e/g' > results1.h
;;   emacs -batch -l test.el -f test2 | sed  's/00\+e/0e/g' > results.h
;;   emacs -batch -l test.el -f test3 | sed  's/00\+e/0e/g' > results_real.h
;;   emacs -batch -l test.el -f test4 | sed  's/00\+e/0e/g' > results2.h
;;   emacs -batch -l test.el -f test5 | sed  's/00\+e/0e/g' > results_zreal.h
;;
;; Note: this takes a long time to run

;; set Calc to use radians mode, turn off symbolic evaluation and use
;; a reasonable precision.

(setq calc-display-working-message t) ;; display short working messages
(setq max-lisp-eval-depth 4000)
(setq max-specpdl-size 2400)
(setq calc-internal-prec 64) ;; probably unnecessarily high, but found
			     ;; a few discrepancies at prec=20
(setq calc-infinite-mode t)
(setq calc-angle-mode 'rad)
(setq calc-float-format '(sci 20))
;;(setq calc-full-float-format '(sci 0))

(setq var-EvalRules (calc-eval "[ sec(x) := 1 / cos(x),
  csc(x) := 1 / sin(x),
  cot(x) := 1 / tan(x),
  sech(x) := 1 / cosh(x),
  csch(x) := 1 / sinh(x),
  coth(x) := 1 / tanh(x),
  arcsec(x) := arccos(1 / x),
  arccsc(x) := arcsin(1 / x),
  arccot(x) := arctan(1 / x),
  arcsech(x) := arccosh(1 / x),
  arccsch(x) := arcsinh(1 / x),
  arccoth(x) := arctanh(1 / x),
  abs2(x) := x * conj(x),
  logabs(x) := log(abs(x)),
  sqrtreal(x) := sqrt(x),
  arcsinreal(x) := arcsin(x),
  arccosreal(x) := arccos(x),
  arccoshreal(x) := arccosh(x),
  arctanhreal(x) := arctanh(x),
  arccscreal(x) := arccsc(x),
  arcsecreal(x) := arcsec(x),
  powreal(x,y) := pow(x,y) ]" 'raw))

;; Convert floating point numbers into exact binary representation

(defun binary (x)
  (if (= x 0)
      x
    (let ((y (/ x (expt 2.0 (+ 1 (logb (abs x)))))))
      (concat (if (>= y 0) "+" "-") (mantissa (abs y)) "e" (format "%d" (+ 1 (logb (abs x))))))))

(defun mantissa (x)
  (let ((y "2#0."))
    (while (> x 0)
      (progn ;(message "x = %g  y = %s" x y)
             (setq x (* 2 x))
             (if (>= x 1) 
                 (progn (setq y (concat y "1"))
                        (setq x (- x 1)))
               (setq y (concat y "0")))))
    y))



;;(binary 9.9999999999999995e-08)


(defun reflections (a b)
  (let ((a (float a)) (b (float b)))
    (list
     (list a b)
     (list a (- b))
     (list (- a) b)
     (list (- a) (- b)))))

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


(defun calcfmt-complex (arg)
   (let* ((x (nth 0 arg))
          (y (nth 1 arg)))
     (format "(%s,%s)" (binary x) (binary y))))

(defun calcfmt-real (arg)
  (format "%s" (binary arg)))

(defun calcfmt (arg)
  (if (listp arg) (calcfmt-complex arg) 
    (calcfmt-real arg)))
      
(defun clean (result)
  (if (string-match "clean(\\(.*\\))" result)
      (setq result (replace-match "\\1" nil nil result)))
  (if (string-match "clean(\\(.*\\), *[0-9]*)" result)
      (setq result (replace-match "\\1" nil nil result)))
  (if (string-match "(\\(.*\\),\\(.*\\))" result)
      (setq result (replace-match "\\1,\\2" nil nil result)))
  (if (string-match "^\\([^,]*\\)$" result)
      (setq result (replace-match "\\1, 0.0" nil nil result)))
  (if (not (or (string-match "(" result) ;; skip any unsimplified results
               (string-match "inf" result)))
      result))

(defun calceval (function arg)
  (let* ((fn (if (string-match "_" function)
                 (replace-match "" nil nil function)
               function))
         (v (concat "0.0 + clean(" fn "(" arg "),60)")))
    ;;(message "v = <%s>" v)
    (clean (calc-eval v))))

(defun cfmt (arg)
  (if (listp arg) (format "ARG(%.20e,%.20e)" (nth  0 arg) (nth  1 arg)) 
    (format "%.20e" arg)))
  
(defun evaltest1 (function arg)
  (let* ((z (calcfmt arg))
         (result (calceval function z)))
    (if result (princ (format "  {FN (%s), %s, RES(%s)},\n" function (cfmt arg) result)))))

(defun evaltest2 (function arg1 arg2)
  (let* ((z (calcfmt arg1))
         (z2 (calcfmt arg2))
         (result (calceval function (concat z "," z2))))
    ;;(message "z=<%s> z2=<%s> result =  <%s>" z z2 result )
    (if result (princ (format "  {FN (%s), %s, %s, RES(%s)},\n" function (cfmt arg1) (cfmt arg2) result)))))

;;(evaltest "sin" "10" "0")

;; loop over all possible combinations of a,b,c

(defun loop1 (a b)
  (let ((b1 b))
    (while b1
      (progn
        (let* ((z (car b1)))
          (evaltest1 a z))
          (setq b1 (cdr b1))
        )
      )
    )
  )

(defun loop2 (a b c)
  (let ((b1 b))
    (while b1
      (progn
        (let* ((z (car b1)) (c1 c))
          (while c1
            (progn 
              (let* ((zp (car c1)))
                (evaltest2 a z zp))
              (setq c1 (cdr c1)))))
          (setq b1 (cdr b1))))))
 
;;(testzz "xx" '((0 1) (2 3)) '((4 5) (6 7)))

;;

(setq pi 3.14159265358979323846264338328);
(setq flteps 1.1920928955078125e-07);
(setq delta (sqrt flteps))

(setq eps (list flteps))
(setq zero (list 0.0))
(setq simple (list 0.5  1.0  2.0))
(setq inf  (list (/ 1.0 flteps)))
(setq finite (append eps simple))
(setq zfinite (append zero eps simple))

(setq realpos (list flteps delta  0.125 0.5 0.75  1.0  2.0  10.0))
(setq real (append (reverse (mapcar '- realpos)) 
                   zero 
                   realpos))

(setq circ (list (- 0 delta)
                 (+ 0 delta)
                 (- (* 0.5 pi) delta)
                 (+ (* 0.5 pi) delta)
                 (- pi delta)
                 (+ pi delta)
                 (- (* 1.5 pi) delta)
                 (+ (* 1.5 pi) delta)
                 (- (* 2 pi) delta)
                 (+ (* 2 pi) delta)
                 (- (* 3 pi) delta)
                 (+ (* 3 pi) delta)))

(setq trig (list (+ (sqrt pi) delta)
                 (+ (log pi) delta)
                 (- (sqrt pi) delta)
                 (- (log pi) delta)))


(setq z2finite (combine (append eps simple) (append zero eps simple)))
(setq z2 (combine (append zero eps simple) (append zero eps simple)))

(setq z0 (combine (append zero eps simple inf)
                  (append zero eps simple inf)))


(setq z1 (append (combine (append eps simple inf)
                          (append zero eps simple inf))
                 (combine (append zero eps simple inf)
                          (append eps simple inf))))


;(setq z2 (append (combine simple real) 
;                 (combine edge (append edge simple) )))

(setq zcirc (append (combine circ zfinite)))
(setq zicirc (append (combine zfinite circ)))

(defun test1 ()
  (loop1 "arg" z0)
  (loop1 "abs" z0)
  (loop1 "abs2" z0)
  (loop1 "logabs" z1)
)

(defun test2 ()
  (loop1 "sqrt" z0)

  (loop1 "log" z1)
  (loop1 "log10" z1)
  (loop1 "exp" zicirc)

  (loop1 "sin" zcirc)
  (loop1 "cos" zcirc)
  (loop1 "tan" zcirc)

  (loop1 "arcsin" z0)
  (loop1 "arccos" z0)
  (loop1 "arctan" z0)

  (loop1 "sinh" zicirc)
  (loop1 "cosh" zicirc)
  (loop1 "tanh" zicirc)

  (loop1 "arcsinh" z0)
  (loop1 "arccosh" z0)
  (loop1 "arctanh" z0)

  (loop1 "csc" zcirc)
  (loop1 "sec" zcirc)
  (loop1 "cot" zcirc)

  (loop1 "arccsc" z0)
  (loop1 "arcsec" z0)
  (loop1 "arccot" z0)

  (loop1 "csch" zicirc)
  (loop1 "sech" zicirc)
  (loop1 "coth" zicirc)

  (loop1 "arccsch" z0)
  (loop1 "arcsech" z0)
  (loop1 "arccoth" z0)
)

(defun test3 ()
  (loop1 "sqrt_real" real)

  (loop1 "arcsin_real" real)
  (loop1 "arccos_real" real)

  (loop1 "arccosh_real" real)
  (loop1 "arctanh_real" real)

  (loop1 "arccsc_real" real)
  (loop1 "arcsec_real" real)
)


(defun test4 ()
  (loop2 "pow" 
         (list '(1 0) '(0 1) '(-1 0) '(0 -1) '(0.5 0.1) '(0.5 -0.1))
         (list '(0 0)  '(1 0) '(0 1) '(-1 0) '(0 -1) '(0.5 0.1) '(0.5 -0.1)))
)

(defun test5 ()
  (loop2 "pow_real" 
         (list '(0 0) '(1 0) '(0 1) '(-1 0) '(0 -1) '(0.5 0.1) '(0.5 -0.1))
         (list 0 1 0.5 2))
)

;;(test1)
;;(test-all)
;;(test3)

;;(evaltestzz "pow" (list flteps flteps) (list -2.0 0.0))
;;(evaltest2 "pow_real" '(1 2) 3)
