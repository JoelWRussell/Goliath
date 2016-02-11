;; gorilla-repl.fileformat = 1

;; **
;;; # Testing JLink
;;; 
;;; You can use this worksheet to test whether JLink is set up correctly on your computer.
;; **

;; @@
(ns warm-swamp
  (:require [gorilla-plot.core :as plot])
  (:import [com.wolfram.jlink]))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@
(def kernel (com.wolfram.jlink.MathLinkFactory/createKernelLink "-linkmode launch -linkname \"/Applications/Mathematica.app/Contents/MacOS/MathKernel\""))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;warm-swamp/kernel</span>","value":"#'warm-swamp/kernel"}
;; <=

;; @@
(.waitForAnswer kernel)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>8</span>","value":"8"}
;; <=

;; @@
(.getString kernel)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-string'>&quot;In[1]:= &quot;</span>","value":"\"In[1]:= \""}
;; <=

;; @@
(.evaluate kernel "1 + 1")
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@
(.waitForAnswer kernel)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>3</span>","value":"3"}
;; <=

;; @@
(.getString kernel)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-string'>&quot;2&quot;</span>","value":"\"2\""}
;; <=

;; @@
(.close kernel)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@

;; @@
