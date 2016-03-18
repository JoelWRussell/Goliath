;; gorilla-repl.fileformat = 1

;; **
;;; # Gorilla REPL
;;; 
;;; Welcome to gorilla :-)
;;; 
;;; Shift + enter evaluates code. Hit alt+g twice in quick succession or click the menu icon (upper-right corner) for more commands ...
;;; 
;;; It's a good habit to run each worksheet in its own namespace: feel free to use the declaration we've provided below if you'd like.
;; **

;; @@
(ns spacial-lake
  (:require [gorilla-plot.core :as plot]))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@
(def indvs [{:a 1 :b 2 :g [1 1 2 2]} {:a 3 :b 4 :g [2 2 1 1]} {:a 2 :b 3 :g [1 2 1 2]}])(map)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;spacial-lake/indvs</span>","value":"#'spacial-lake/indvs"}
;; <=

;; @@
(def gen (map :g indvs))
gen
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"}],"value":"[1 1 2 2]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"}],"value":"[2 2 1 1]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"}],"value":"[1 2 1 2]"}],"value":"([1 1 2 2] [2 2 1 1] [1 2 1 2])"}
;; <=

;; @@
(def compl (map count gen))
compl
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"},{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"},{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"}],"value":"(4 4 4)"}
;; <=

;; @@
(def error (map #(reduce + %) gen))
error
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>6</span>","value":"6"},{"type":"html","content":"<span class='clj-long'>6</span>","value":"6"},{"type":"html","content":"<span class='clj-long'>6</span>","value":"6"}],"value":"(6 6 6)"}
;; <=

;; @@
(apply hash-map [:a 1 :b 2])
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-map'>{</span>","close":"<span class='clj-map'>}</span>","separator":", ","items":[{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:b</span>","value":":b"},{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"}],"value":"[:b 2]"},{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:a</span>","value":":a"},{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"}],"value":"[:a 1]"}],"value":"{:b 2, :a 1}"}
;; <=

;; @@
(map #(apply hash-map [:compl %1 :error %2]) compl error)
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-map'>{</span>","close":"<span class='clj-map'>}</span>","separator":", ","items":[{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:compl</span>","value":":compl"},{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"}],"value":"[:compl 4]"},{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:error</span>","value":":error"},{"type":"html","content":"<span class='clj-long'>6</span>","value":"6"}],"value":"[:error 6]"}],"value":"{:compl 4, :error 6}"},{"type":"list-like","open":"<span class='clj-map'>{</span>","close":"<span class='clj-map'>}</span>","separator":", ","items":[{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:compl</span>","value":":compl"},{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"}],"value":"[:compl 4]"},{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:error</span>","value":":error"},{"type":"html","content":"<span class='clj-long'>6</span>","value":"6"}],"value":"[:error 6]"}],"value":"{:compl 4, :error 6}"},{"type":"list-like","open":"<span class='clj-map'>{</span>","close":"<span class='clj-map'>}</span>","separator":", ","items":[{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:compl</span>","value":":compl"},{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"}],"value":"[:compl 4]"},{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:error</span>","value":":error"},{"type":"html","content":"<span class='clj-long'>6</span>","value":"6"}],"value":"[:error 6]"}],"value":"{:compl 4, :error 6}"}],"value":"({:compl 4, :error 6} {:compl 4, :error 6} {:compl 4, :error 6})"}
;; <=

;; @@
(defn fff [indvs]
  
  (let[gen (map :genome indvs)
      compl (map count gen)
      error (map #(reduce + %) gen)
  	  values (map #(apply hash-map [:complexity %1 :error %2]) compl error)]
   (map #(merge %1 %2) indvs values)
  ))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;spacial-lake/fff</span>","value":"#'spacial-lake/fff"}
;; <=

;; @@
(fff indvs)
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-lazy-seq'>(</span>","close":"<span class='clj-lazy-seq'>)</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-map'>{</span>","close":"<span class='clj-map'>}</span>","separator":", ","items":[{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:complexity</span>","value":":complexity"},{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"}],"value":"[:complexity 4]"},{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:g</span>","value":":g"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"}],"value":"[1 1 2 2]"}],"value":"[:g [1 1 2 2]]"},{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:b</span>","value":":b"},{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"}],"value":"[:b 2]"},{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:error</span>","value":":error"},{"type":"html","content":"<span class='clj-long'>6</span>","value":"6"}],"value":"[:error 6]"},{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:a</span>","value":":a"},{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"}],"value":"[:a 1]"}],"value":"{:complexity 4, :g [1 1 2 2], :b 2, :error 6, :a 1}"},{"type":"list-like","open":"<span class='clj-map'>{</span>","close":"<span class='clj-map'>}</span>","separator":", ","items":[{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:complexity</span>","value":":complexity"},{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"}],"value":"[:complexity 4]"},{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:g</span>","value":":g"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"}],"value":"[2 2 1 1]"}],"value":"[:g [2 2 1 1]]"},{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:b</span>","value":":b"},{"type":"html","content":"<span class='clj-long'>4</span>","value":"4"}],"value":"[:b 4]"},{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:error</span>","value":":error"},{"type":"html","content":"<span class='clj-long'>6</span>","value":"6"}],"value":"[:error 6]"},{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:a</span>","value":":a"},{"type":"html","content":"<span class='clj-long'>3</span>","value":"3"}],"value":"[:a 3]"}],"value":"{:complexity 4, :g [2 2 1 1], :b 4, :error 6, :a 3}"},{"type":"list-like","open":"<span class='clj-map'>{</span>","close":"<span class='clj-map'>}</span>","separator":", ","items":[{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:complexity</span>","value":":complexity"},{"type":"html","content":"<span class='clj-unkown'>4</span>","value":"4"}],"value":"[:complexity 4]"},{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:g</span>","value":":g"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"}],"value":"[1 2 1 2]"}],"value":"[:g [1 2 1 2]]"},{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:b</span>","value":":b"},{"type":"html","content":"<span class='clj-long'>3</span>","value":"3"}],"value":"[:b 3]"},{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:error</span>","value":":error"},{"type":"html","content":"<span class='clj-long'>6</span>","value":"6"}],"value":"[:error 6]"},{"type":"list-like","open":"","close":"","separator":" ","items":[{"type":"html","content":"<span class='clj-keyword'>:a</span>","value":":a"},{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"}],"value":"[:a 2]"}],"value":"{:complexity 4, :g [1 2 1 2], :b 3, :error 6, :a 2}"}],"value":"({:complexity 4, :g [1 1 2 2], :b 2, :error 6, :a 1} {:complexity 4, :g [2 2 1 1], :b 4, :error 6, :a 3} {:complexity 4, :g [1 2 1 2], :b 3, :error 6, :a 2})"}
;; <=

;; @@

;; @@
