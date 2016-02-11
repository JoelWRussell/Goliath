;; gorilla-repl.fileformat = 1

;; **
;;; # Link To Polynomialiser
;;; This Gorilla worksheet shows how a polynomial can be scored
;;; 
;; **

;; @@
(ns lagrangianPolynomialiser
  (:require [gorilla-plot.core :as plot]))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; **
;;; First load the jar lagrangianscore - this is fetched in via the Project.clj
;; **

;; @@
(import 'goliath.mathlink.LagrangianScore)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-class'>goliath.mathlink.LagrangianScore</span>","value":"goliath.mathlink.LagrangianScore"}
;; <=

;; **
;;; InitFunctions starts the Mathematica kernel and also loads the pendulum data; it also loads all of the score functions into the kernel.
;; **

;; @@
(LagrangianScore/InitFunctions)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; **
;;; Demonstrate with 2 polnomials poly1 and poly2
;; **

;; **
;;; 
;; **

;; **
;;; 
;; **

;; @@
(def poly1 [1,1,1,1, 2,2,2,2])
(def poly2 [1,1,1,1, 1,2,1,2, 2,1,1,2])

;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;lagrangianPolynomialiser/poly2</span>","value":"#'lagrangianPolynomialiser/poly2"}
;; <=

;; **
;;; Now score each one; the result is the form [score, c0,c1,c2,...]
;; **

;; @@
(into [] (LagrangianScore/GetScore (into-array Integer/TYPE poly1)))
(into [] (LagrangianScore/GetScore (into-array Integer/TYPE poly2)))
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>18.429935563043156</span>","value":"18.429935563043156"},{"type":"html","content":"<span class='clj-double'>1.0</span>","value":"1.0"},{"type":"html","content":"<span class='clj-double'>0.052252465404158316</span>","value":"0.052252465404158316"},{"type":"html","content":"<span class='clj-double'>0.05033719645452601</span>","value":"0.05033719645452601"}],"value":"[18.429935563043156 1.0 0.052252465404158316 0.05033719645452601]"}
;; <=

;; **
;;; Now close the Mathematica kernel
;; **

;; @@
(LagrangianScore/Shutdown);
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@

;; @@
