;; gorilla-repl.fileformat = 1

;; **
;;; #HYBRID WITH VARIABLE DEGREES OF FREEDOM
;;; 
;;; Join up the GP search with the Mathematica search. You can also enter the number of independent df. Assuming that each column in the data set is a different coordinate.
;; **

;; @@
(ns goliath
  (:require [gorilla-plot.core :as plot]
            [clojure.repl :as repl]
            [clojure.pprint :as pprint]
            [clojure.walk :as walk]
            [clojure.zip :as zip]
            [gorilla-repl.table :as table]
            [gorilla-repl.html :as html]
    		[darwin.core :as darwin]
            [darwin.evolution.metrics :as metrics]
            [algebolic.expression.core :as expression]
            [algebolic.expression.tree :as tree]
            [algebolic.expression.genetics :as genetics]
            [algebolic.expression.score :as score]
            [algebolic.expression.render :as spea-render]
            [algebolic.expression.interpreter :as interpreter]
            [darwin.evolution.core :as evolution]
            [darwin.evolution.metrics :as metrics]
            [darwin.evolution.reproduction :as reproduction]
            [darwin.evolution.scoring :as scoring]
            [darwin.evolution.selection :as selection]
            [darwin.evolution.transform :as transform]
            [darwin.evolution.pareto :as pareto]
            [darwin.algorithms.spea2 :as spea2]
            [criterium.core :as criterium])
  (:import [com.wolfram.jlink]
      [goliath.mathlink])
  
  
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; **
;;; SET THE PATH TO THE MATHEMATICA KERNEL
;; **

;; @@
(def mathKernelSz "c:/program files/wolfram research/mathematica/10.0/mathkernel.exe")
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/mathKernelSz</span>","value":"#'goliath/mathKernelSz"}
;; <=

;; **
;;; CHOOSE THE DATA - NEEDS TO BE IN RESOURCES/
;; **

;; @@
(def experimentalDataSz "mma_coupled_sho.csv")
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/experimentalDataSz</span>","value":"#'goliath/experimentalDataSz"}
;; <=

;; **
;;; CHOOSE THE TIME INTERVAL
;; **

;; @@
(def dt 0.1)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/dt</span>","value":"#'goliath/dt"}
;; <=

;; **
;;; CHOOSE THE NUMBER OF DEGREES OF FREEDOM (in this case 2 for the double prendulum)
;; **

;; @@
(def df 2)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/df</span>","value":"#'goliath/df"}
;; <=

;; @@
(def prob_inheritance 0.75)   ;probability that terms are directly inherited from parent to child.


(def mutate_pref1 0.9); (1 - mutate_pref1) = prob of gene-delete     
(def mutate_pref2 0.8); (mutate_pref1 - mutate_pref2) = prob. of gene-tweak
(def mutate_pref3 0.4);(mutate_pref2 - mutate_pref3) = prob. of gene-replace,
					  ; mutate_pref3 = prob. of gene-add
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/mutate_pref3</span>","value":"#'goliath/mutate_pref3"}
;; <=

;; @@
(def power_sum 10)    ;the maximum sum of powers present in each term within a polynomial.
(def poly_length 10)  ;the maximum number of terms of a polynomial. 
(def power 8) 		  ;the maximum power a variable can be raised to.
(def term_length (* 2 df))   ;the number of variables within each term - 2*degree_freedom
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/term_length</span>","value":"#'goliath/term_length"}
;; <=

;; @@
(defn new-poly-term 
  
  "creates a random new poly term satisfying given constraints
   powermax - the maximum power of any given variable, 
   spmax - the max sum of the powers and length is the length of the gene to be generated. 
   spmax - checked and the function re-called if the gene is invalid"
  
  [max_pow spmax length]
  
  (let [new-gene (vec (repeatedly length #(rand-int (+ 1 max_pow)))) sp (apply + new-gene)]
    
    (if (or (< spmax sp) (= new-gene (vec (replicate length 0)) )) 
      
      (assoc (vec (repeat length 0)) (rand-int length) 1)
    
      ;(new-poly-term max_pow spmax length)
      ;[0 1];;cant have this or else it will crash the mathematica
      new-gene
      )
  ))
  
  
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/new-poly-term</span>","value":"#'goliath/new-poly-term"}
;; <=

;; @@
(new-poly-term 10 3 3)
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-long'>0</span>","value":"0"},{"type":"html","content":"<span class='clj-long'>0</span>","value":"0"}],"value":"[1 0 0]"}
;; <=

;; @@
(defn create-genotype
  
  "creates a new individual of length between 1 and max_length decided at random.
  terms within an individual are of legnth n.
  polynomial terms are constrained to maximum power p_max"
  
  [max_length max_pow spmax term_length]
  
  (let [geno (vec (repeatedly (+ 1 (rand-int (- max_length 1))) #(new-poly-term max_pow spmax term_length)))]
    
      (vec (distinct geno))
     )    
    )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/create-genotype</span>","value":"#'goliath/create-genotype"}
;; <=

;; @@
(defn distribute-parent   
  "Distributes the alleses of a parent to children with seeded probability"
 ([parent-1 prob]
  (let [{c1 true c2 false} (group-by (fn [x] (<= prob (rand))) parent-1)]    
[c1 c2]  
   ))
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/distribute-parent</span>","value":"#'goliath/distribute-parent"}
;; <=

;; **
;;; 
;; **

;; **
;;; ##Cross-over and mutate functions
;; **

;; **
;;; cross-over, gene-mutate, gene-replace and gene-tweak
;; **

;; @@
(defn cross-over
  
  "takes two individuals and crosses over the polynomial from one parent to another under a weighted selection bias" 
  
  [parent-1 parent-2] 
  
  (let [offspring (mapv vec (mapv set (mapv concat (distribute-parent parent-1 prob_inheritance) (reverse (distribute-parent parent-2 prob_inheritance)))))]
    
    (if (or (zero? (count (offspring 0))) (zero? (count (offspring 1))))
      	
      [parent-1 parent-2]
      [(offspring 0)(offspring 1)]
      )
    ))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/cross-over</span>","value":"#'goliath/cross-over"}
;; <=

;; @@
(defn gene-replace
  
  "Replaces a gene from an individual (with a unique gene). In a given gene (term in the polynomial) no variable can have a higher power than powermax and the total powers of all variables cannot exceed spmax (sum power max)."
  
  [indv powermax spmax n]
  
(let [length (count (first indv)) new-gene (new-poly-term powermax spmax length) new-geno (assoc indv (rand-int (count indv)) new-gene) sp (apply + new-gene)]
  
  (if (> n 10) indv  ;;
  (if (= (count (set new-geno)) (count new-geno))
    new-geno
    (gene-replace indv powermax spmax (inc n))
 
 ))))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/gene-replace</span>","value":"#'goliath/gene-replace"}
;; <=

;; @@
(defn gene-add
  
  "Adds a new gene to a GA individual. The gene will be unique addition to the individual. In a given gene (term in the polynomial) no variable can have a higher power than powermax and the total powers of all variables cannot exceed spmax (sum power max)."
  
  [indiv powermax spmax n]
   
  (let [length (count (first indiv)) new-gene (new-poly-term powermax spmax length) new-geno (vec (conj (set indiv) new-gene)) sp (apply + new-gene)]
   
      (if (> n 10) 
        indiv (if (= (count new-geno) (count indiv))
      (gene-add indiv powermax spmax (inc n))
     new-geno
      ))))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/gene-add</span>","value":"#'goliath/gene-add"}
;; <=

;; @@
(defn gene-tweak
  [indv]
  (let [n (rand-int (count indv)) gene (nth indv n) new-gene (assoc gene (rand-int (count gene)) (rand-int power))]
    (assoc indv n new-gene)
    ))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/gene-tweak</span>","value":"#'goliath/gene-tweak"}
;; <=

;; @@
(defn gene-delete
  [indv]
  
  (let [new-indv (subvec indv 1)]
    (if (= new-indv []) indv new-indv))
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/gene-delete</span>","value":"#'goliath/gene-delete"}
;; <=

;; @@
(gene-delete [[10 0 0 0] [5 4 3 2]])
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>5</span>","value":"5"},{"type":"html","content":"<span class='clj-long'>4</span>","value":"4"},{"type":"html","content":"<span class='clj-long'>3</span>","value":"3"},{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"}],"value":"[5 4 3 2]"}],"value":"[[5 4 3 2]]"}
;; <=

;; @@
(defn mutate
  "takes two genotypes and returns a mutated genotype"
  [indv]
  (let [rn (rand)]
    
    (if (>= rn mutate_pref1) 
      (gene-delete indv) 
    (if (>= rn mutate_pref2)
        (gene-tweak indv)
    (if (>= rn mutate_pref3) 
        (gene-replace indv power power_sum 1)
          (gene-add indv power power_sum 1))
      )
     ) 
    )
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/mutate</span>","value":"#'goliath/mutate"}
;; <=

;; **
;;; #Initial populations and generation-configuration
;; **

;; **
;;; plus a little kernal opening and closing
;; **

;; @@
(defn random-initial-population 
  
   "creates an innitial population of polynimials
  size - number of individuals in the population"
   
  [size max_length max_power max_sum_powers num_terms]
 
  (repeatedly size #(create-genotype max_length max_power max_sum_powers num_terms)) 
       
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/random-initial-population</span>","value":"#'goliath/random-initial-population"}
;; <=

;; @@
(def initial-zeitgeist (evolution/make-zeitgeist (random-initial-population 100 poly_length power power_sum term_length)))

;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/initial-zeitgeist</span>","value":"#'goliath/initial-zeitgeist"}
;; <=

;; @@
(goliath.mathlink.LagrangianScore/InitFunctions
  mathKernelSz
  "resources/"
  experimentalDataSz
  dt
  df
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@
(defn score
  [indv] 
	(first
      (into [] (goliath.mathlink.LagrangianScore/GetScore (into-array Integer/TYPE (vec (flatten indv))), df))
      )
  
  
  )
  

;;(score (create-genotype poly_length power power_sum term_length) )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/score</span>","value":"#'goliath/score"}
;; <=

;; @@
(def generation-config			      
  (let [size-limit 120                   
        min-size 20
        ea-config (spea2/spea2-config
                    {:goals [:error :complexity]
                     :archive-size 50
                     :binary-ops [{:op cross-over :repeat 10}]
                     :unary-ops [{:op mutate :repeat 30}]})
        score-functions {:complexity (fn [x] (count x))
                         :error (fn [e] (score e))}]
    {:ea-config              ea-config
     :score-functions        score-functions
     :reporting-function     (fn [z] (print ".") (flush))}))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/generation-config</span>","value":"#'goliath/generation-config"}
;; <=

;; @@
(time (def result (evolution/run-evolution generation-config initial-zeitgeist (fn [zg gc] (>= (:age zg) 10)))))
;; @@
;; ->
;;; ..........&quot;Elapsed time: 17378.824395 msecs&quot;
;;; 
;; <-
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/result</span>","value":"#'goliath/result"}
;; <=

;; @@
(goliath.mathlink.LagrangianScore/Shutdown)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@
(mapv #(println (:genotype %)) (sort-by :error (:elite result)))
;; @@
;; ->
;;; [[4 2 0 4] [1 0 0 0] [0 1 0 0] [0 0 2 7] [1 1 1 0] [2 2 2 1] [1 0 8 0]]
;;; [[1 0 0 0] [0 1 0 0] [0 0 2 7] [1 1 1 0] [2 3 2 1] [1 0 8 0]]
;;; [[1 0 0 0] [0 1 0 0] [0 0 2 7] [1 1 1 0] [1 0 8 0]]
;;; [[1 0 0 0] [0 0 2 7] [1 1 1 0] [1 0 8 0]]
;;; [[1 0 0 0] [2 2 2 2] [1 1 1 0]]
;;; [[1 0 0 0] [1 1 1 0]]
;;; [[1 0 0 0] [1 1 1 0]]
;;; [[1 0 0 0] [1 1 1 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[0 1 0 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; 
;; <-
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}],"value":"[nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil]"}
;; <=

;; **
;;; #Pareto Analysis
;;; 
;; **

;; **
;;; s
;; **

;; @@
(plot/list-plot (:min (:error @metrics/metrics)) :joined true)
;; @@
;; =>
;;; {"type":"vega","content":{"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"1930277d-a27e-40fa-8843-8e6b32f844f4","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":{"data":"1930277d-a27e-40fa-8843-8e6b32f844f4","field":"data.y"}}],"marks":[{"type":"line","from":{"data":"1930277d-a27e-40fa-8843-8e6b32f844f4"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"#FF29D2"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}}],"data":[{"name":"1930277d-a27e-40fa-8843-8e6b32f844f4","values":[{"x":0,"y":-0.590430086133494},{"x":1,"y":-0.6164908867191314},{"x":2,"y":-0.6318432531955267},{"x":3,"y":-0.7200761681601696},{"x":4,"y":-0.7423913782757169},{"x":5,"y":-0.913357253037232},{"x":6,"y":-0.913357253037232},{"x":7,"y":-0.913357253037232},{"x":8,"y":-0.913357253037232},{"x":9,"y":-0.913357253037232}]}],"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50}},"value":"#gorilla_repl.vega.VegaView{:content {:axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"1930277d-a27e-40fa-8843-8e6b32f844f4\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain {:data \"1930277d-a27e-40fa-8843-8e6b32f844f4\", :field \"data.y\"}}], :marks [{:type \"line\", :from {:data \"1930277d-a27e-40fa-8843-8e6b32f844f4\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"#FF29D2\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}}], :data [{:name \"1930277d-a27e-40fa-8843-8e6b32f844f4\", :values ({:x 0, :y -0.590430086133494} {:x 1, :y -0.6164908867191314} {:x 2, :y -0.6318432531955267} {:x 3, :y -0.7200761681601696} {:x 4, :y -0.7423913782757169} {:x 5, :y -0.913357253037232} {:x 6, :y -0.913357253037232} {:x 7, :y -0.913357253037232} {:x 8, :y -0.913357253037232} {:x 9, :y -0.913357253037232})}], :width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}}}"}
;; <=

;; @@
(defn pareto-plot-population
  [[k1 k2] result]
  (let [pareto-front (pareto/non-dominated-individuals [k1 k2] (:elite result))
        coord-extract (fn [i] [(k1 i) (k2 i)])]
    (plot/compose
      (plot/list-plot (map coord-extract (:elite result)) :colour "red")
      (plot/list-plot (map coord-extract (:rabble result)) :colour "blue")
          (plot/list-plot (map coord-extract pareto-front) :colour "#ff29d2")
  )))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/pareto-plot-population</span>","value":"#'goliath/pareto-plot-population"}
;; <=

;; @@
(pareto-plot-population [:error :complexity] result)
;; @@
;; =>
;;; {"type":"vega","content":{"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50},"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"fbec7a07-1e35-4738-8f6f-4227b05f895e","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":{"data":"fbec7a07-1e35-4738-8f6f-4227b05f895e","field":"data.y"}}],"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"data":[{"name":"fbec7a07-1e35-4738-8f6f-4227b05f895e","values":[{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":-0.5884356859669099,"y":2},{"x":-0.5884356859669099,"y":2},{"x":-0.5884356859669099,"y":2},{"x":-0.6629791795518623,"y":5},{"x":-0.6636021818967084,"y":6},{"x":-0.6620521442016984,"y":4},{"x":-0.6203403511070944,"y":3},{"x":-0.913357253037232,"y":7}]},{"name":"74c2adc1-e5b4-4e65-8226-a2b1394d916c","values":[{"x":-1.2667705298363256E-8,"y":2},{"x":14.311479305036924,"y":1},{"x":16.212236546842988,"y":1},{"x":-0.663010216898583,"y":6},{"x":16.212236546842988,"y":1},{"x":0.00986814613767775,"y":2},{"x":-0.6211468434957376,"y":4},{"x":1.5832220019188459,"y":1},{"x":-2.3889364493379955E-6,"y":2},{"x":-0.7059460269039661,"y":7},{"x":1.5832220019188459,"y":1},{"x":0.0,"y":2},{"x":1.5832220019188459,"y":1},{"x":-1.2707238978988279E-6,"y":2},{"x":-0.5884356859669099,"y":3},{"x":9.07872731939477,"y":1},{"x":9.07872731939477,"y":2},{"x":-0.58860558252265,"y":3},{"x":16.212236546842988,"y":1},{"x":-0.0028376297044913333,"y":5},{"x":-5.645022542681962E-5,"y":2},{"x":0.00986814613767775,"y":2},{"x":0.00986814613767775,"y":2},{"x":0.00986814613767775,"y":2},{"x":9.07872731939477,"y":2},{"x":0.0,"y":2},{"x":-0.6211468434957376,"y":4},{"x":16.212236546842988,"y":1},{"x":0.00986814613767775,"y":2},{"x":1.5832220019188459,"y":1},{"x":-0.6629791795518623,"y":5},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":9.07872731939477,"y":1},{"x":0.0,"y":2},{"x":-0.5884356859669099,"y":2},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":-0.6629791795518623,"y":5},{"x":-0.0027840623537830488,"y":3},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":-0.7412925940852366,"y":5},{"x":36.554339722388434,"y":2}]},{"name":"39d33455-b647-463e-8414-28e7a0ffd416","values":[{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":1.5832220019188459,"y":1},{"x":-0.5884356859669099,"y":2},{"x":-0.5884356859669099,"y":2},{"x":-0.5884356859669099,"y":2},{"x":-0.6629791795518623,"y":5},{"x":-0.6636021818967084,"y":6},{"x":-0.6620521442016984,"y":4},{"x":-0.6203403511070944,"y":3},{"x":-0.913357253037232,"y":7}]}],"marks":[{"type":"symbol","from":{"data":"fbec7a07-1e35-4738-8f6f-4227b05f895e"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"fill":{"value":"red"},"fillOpacity":{"value":1}},"update":{"shape":"circle","size":{"value":70},"stroke":{"value":"transparent"}},"hover":{"size":{"value":210},"stroke":{"value":"white"}}}},{"type":"symbol","from":{"data":"74c2adc1-e5b4-4e65-8226-a2b1394d916c"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"fill":{"value":"blue"},"fillOpacity":{"value":1}},"update":{"shape":"circle","size":{"value":70},"stroke":{"value":"transparent"}},"hover":{"size":{"value":210},"stroke":{"value":"white"}}}},{"type":"symbol","from":{"data":"39d33455-b647-463e-8414-28e7a0ffd416"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"fill":{"value":"#ff29d2"},"fillOpacity":{"value":1}},"update":{"shape":"circle","size":{"value":70},"stroke":{"value":"transparent"}},"hover":{"size":{"value":210},"stroke":{"value":"white"}}}}]},"value":"#gorilla_repl.vega.VegaView{:content {:width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}, :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"fbec7a07-1e35-4738-8f6f-4227b05f895e\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain {:data \"fbec7a07-1e35-4738-8f6f-4227b05f895e\", :field \"data.y\"}}], :axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :data ({:name \"fbec7a07-1e35-4738-8f6f-4227b05f895e\", :values ({:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x -0.5884356859669099, :y 2} {:x -0.5884356859669099, :y 2} {:x -0.5884356859669099, :y 2} {:x -0.6629791795518623, :y 5} {:x -0.6636021818967084, :y 6} {:x -0.6620521442016984, :y 4} {:x -0.6203403511070944, :y 3} {:x -0.913357253037232, :y 7})} {:name \"74c2adc1-e5b4-4e65-8226-a2b1394d916c\", :values ({:x -1.2667705298363256E-8, :y 2} {:x 14.311479305036924, :y 1} {:x 16.212236546842988, :y 1} {:x -0.663010216898583, :y 6} {:x 16.212236546842988, :y 1} {:x 0.00986814613767775, :y 2} {:x -0.6211468434957376, :y 4} {:x 1.5832220019188459, :y 1} {:x -2.3889364493379955E-6, :y 2} {:x -0.7059460269039661, :y 7} {:x 1.5832220019188459, :y 1} {:x 0.0, :y 2} {:x 1.5832220019188459, :y 1} {:x -1.2707238978988279E-6, :y 2} {:x -0.5884356859669099, :y 3} {:x 9.07872731939477, :y 1} {:x 9.07872731939477, :y 2} {:x -0.58860558252265, :y 3} {:x 16.212236546842988, :y 1} {:x -0.0028376297044913333, :y 5} {:x -5.645022542681962E-5, :y 2} {:x 0.00986814613767775, :y 2} {:x 0.00986814613767775, :y 2} {:x 0.00986814613767775, :y 2} {:x 9.07872731939477, :y 2} {:x 0.0, :y 2} {:x -0.6211468434957376, :y 4} {:x 16.212236546842988, :y 1} {:x 0.00986814613767775, :y 2} {:x 1.5832220019188459, :y 1} {:x -0.6629791795518623, :y 5} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 9.07872731939477, :y 1} {:x 0.0, :y 2} {:x -0.5884356859669099, :y 2} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x -0.6629791795518623, :y 5} {:x -0.0027840623537830488, :y 3} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x -0.7412925940852366, :y 5} {:x 36.554339722388434, :y 2})} {:name \"39d33455-b647-463e-8414-28e7a0ffd416\", :values ({:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x 1.5832220019188459, :y 1} {:x -0.5884356859669099, :y 2} {:x -0.5884356859669099, :y 2} {:x -0.5884356859669099, :y 2} {:x -0.6629791795518623, :y 5} {:x -0.6636021818967084, :y 6} {:x -0.6620521442016984, :y 4} {:x -0.6203403511070944, :y 3} {:x -0.913357253037232, :y 7})}), :marks ({:type \"symbol\", :from {:data \"fbec7a07-1e35-4738-8f6f-4227b05f895e\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :fill {:value \"red\"}, :fillOpacity {:value 1}}, :update {:shape \"circle\", :size {:value 70}, :stroke {:value \"transparent\"}}, :hover {:size {:value 210}, :stroke {:value \"white\"}}}} {:type \"symbol\", :from {:data \"74c2adc1-e5b4-4e65-8226-a2b1394d916c\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :fill {:value \"blue\"}, :fillOpacity {:value 1}}, :update {:shape \"circle\", :size {:value 70}, :stroke {:value \"transparent\"}}, :hover {:size {:value 210}, :stroke {:value \"white\"}}}} {:type \"symbol\", :from {:data \"39d33455-b647-463e-8414-28e7a0ffd416\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :fill {:value \"#ff29d2\"}, :fillOpacity {:value 1}}, :update {:shape \"circle\", :size {:value 70}, :stroke {:value \"transparent\"}}, :hover {:size {:value 210}, :stroke {:value \"white\"}}}})}}"}
;; <=

;; **
;;; Timings of the various steps in each generation.
;;; 
;;; - "Red" - selection time
;;; - "Green" Reproduction time 
;;; - "Blue"  Scoring time
;; **

;; @@
(plot/compose
  (plot/list-plot (:time @metrics/metrics) :colour "yellow" :joined true
                  :plot-range [:all [0 (apply max (:time @metrics/metrics))]])
  (plot/list-plot (:selection-time @metrics/metrics) :colour "red" :joined true)
  (plot/list-plot (:reproduction-time @metrics/metrics) :colour "green" :joined true)
  (plot/list-plot (:scoring-time @metrics/metrics) :colour "blue" :joined true))
;; @@
;; =>
;;; {"type":"vega","content":{"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50},"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"50111b8e-ed9e-496f-aced-d45a43d2f4d3","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":[0,36842]}],"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"data":[{"name":"50111b8e-ed9e-496f-aced-d45a43d2f4d3","values":[{"x":0,"y":10030},{"x":1,"y":12567},{"x":2,"y":5347},{"x":3,"y":7512},{"x":4,"y":11195},{"x":5,"y":15906},{"x":6,"y":11554},{"x":7,"y":15888},{"x":8,"y":18141},{"x":9,"y":14429},{"x":10,"y":19244},{"x":11,"y":18015},{"x":12,"y":19820},{"x":13,"y":32980},{"x":14,"y":29509},{"x":15,"y":20270},{"x":16,"y":28139},{"x":17,"y":25937},{"x":18,"y":20463},{"x":19,"y":30241},{"x":20,"y":27274},{"x":21,"y":28576},{"x":22,"y":32776},{"x":23,"y":30361},{"x":24,"y":21013},{"x":25,"y":22067},{"x":26,"y":19214},{"x":27,"y":28385},{"x":28,"y":23763},{"x":29,"y":30851},{"x":30,"y":12249},{"x":31,"y":9257},{"x":32,"y":11395},{"x":33,"y":9844},{"x":34,"y":16059},{"x":35,"y":18623},{"x":36,"y":16781},{"x":37,"y":26557},{"x":38,"y":20008},{"x":39,"y":20427},{"x":40,"y":21009},{"x":41,"y":26002},{"x":42,"y":10923},{"x":43,"y":8629},{"x":44,"y":10014},{"x":45,"y":13493},{"x":46,"y":29692},{"x":47,"y":35534},{"x":48,"y":31156},{"x":49,"y":21538},{"x":50,"y":19802},{"x":51,"y":15251},{"x":52,"y":9091},{"x":53,"y":7819},{"x":54,"y":8359},{"x":55,"y":8648},{"x":56,"y":8863},{"x":57,"y":8637},{"x":58,"y":8484},{"x":59,"y":11055},{"x":60,"y":11655},{"x":61,"y":14464},{"x":62,"y":13165},{"x":63,"y":15865},{"x":64,"y":17436},{"x":65,"y":17009},{"x":66,"y":22685},{"x":67,"y":24846},{"x":68,"y":18818},{"x":69,"y":30613},{"x":70,"y":21982},{"x":71,"y":18029},{"x":72,"y":16339},{"x":73,"y":26569},{"x":74,"y":23324},{"x":75,"y":33448},{"x":76,"y":36842},{"x":77,"y":30567},{"x":78,"y":31721},{"x":79,"y":15527},{"x":80,"y":15795},{"x":81,"y":21516},{"x":82,"y":23033},{"x":83,"y":19970},{"x":84,"y":22183},{"x":85,"y":18054},{"x":86,"y":26898},{"x":87,"y":26488},{"x":88,"y":16755},{"x":89,"y":18360},{"x":90,"y":20701},{"x":91,"y":17628},{"x":92,"y":22793},{"x":93,"y":21204},{"x":94,"y":25189},{"x":95,"y":23571},{"x":96,"y":23032},{"x":97,"y":20836},{"x":98,"y":21992},{"x":99,"y":26615}]},{"name":"66bce705-a4ba-43ed-bc8a-0d9f1aee5652","values":[{"x":0,"y":64},{"x":1,"y":48},{"x":2,"y":37},{"x":3,"y":26},{"x":4,"y":27},{"x":5,"y":28},{"x":6,"y":26},{"x":7,"y":22},{"x":8,"y":23},{"x":9,"y":20},{"x":10,"y":25},{"x":11,"y":18},{"x":12,"y":19},{"x":13,"y":20},{"x":14,"y":21},{"x":15,"y":20},{"x":16,"y":22},{"x":17,"y":22},{"x":18,"y":18},{"x":19,"y":19},{"x":20,"y":17},{"x":21,"y":16},{"x":22,"y":15},{"x":23,"y":19},{"x":24,"y":24},{"x":25,"y":21},{"x":26,"y":26},{"x":27,"y":28},{"x":28,"y":16},{"x":29,"y":16},{"x":30,"y":24},{"x":31,"y":50},{"x":32,"y":40},{"x":33,"y":42},{"x":34,"y":18},{"x":35,"y":20},{"x":36,"y":19},{"x":37,"y":18},{"x":38,"y":20},{"x":39,"y":23},{"x":40,"y":45},{"x":41,"y":31},{"x":42,"y":19},{"x":43,"y":25},{"x":44,"y":39},{"x":45,"y":37},{"x":46,"y":34},{"x":47,"y":19},{"x":48,"y":18},{"x":49,"y":17},{"x":50,"y":19},{"x":51,"y":24},{"x":52,"y":24},{"x":53,"y":39},{"x":54,"y":41},{"x":55,"y":39},{"x":56,"y":50},{"x":57,"y":27},{"x":58,"y":25},{"x":59,"y":30},{"x":60,"y":29},{"x":61,"y":25},{"x":62,"y":30},{"x":63,"y":28},{"x":64,"y":28},{"x":65,"y":24},{"x":66,"y":33},{"x":67,"y":28},{"x":68,"y":15},{"x":69,"y":30},{"x":70,"y":19},{"x":71,"y":17},{"x":72,"y":25},{"x":73,"y":23},{"x":74,"y":21},{"x":75,"y":24},{"x":76,"y":22},{"x":77,"y":23},{"x":78,"y":28},{"x":79,"y":17},{"x":80,"y":14},{"x":81,"y":28},{"x":82,"y":22},{"x":83,"y":20},{"x":84,"y":17},{"x":85,"y":22},{"x":86,"y":25},{"x":87,"y":29},{"x":88,"y":19},{"x":89,"y":25},{"x":90,"y":26},{"x":91,"y":27},{"x":92,"y":28},{"x":93,"y":19},{"x":94,"y":22},{"x":95,"y":21},{"x":96,"y":19},{"x":97,"y":19},{"x":98,"y":24},{"x":99,"y":20}]},{"name":"d15a96e5-b6a5-4601-9760-399ad56d5c56","values":[{"x":0,"y":116},{"x":1,"y":5},{"x":2,"y":4},{"x":3,"y":2},{"x":4,"y":3},{"x":5,"y":4},{"x":6,"y":2},{"x":7,"y":3},{"x":8,"y":4},{"x":9,"y":2},{"x":10,"y":2},{"x":11,"y":2},{"x":12,"y":2},{"x":13,"y":2},{"x":14,"y":2},{"x":15,"y":2},{"x":16,"y":1},{"x":17,"y":2},{"x":18,"y":3},{"x":19,"y":2},{"x":20,"y":3},{"x":21,"y":2},{"x":22,"y":2},{"x":23,"y":2},{"x":24,"y":2},{"x":25,"y":2},{"x":26,"y":3},{"x":27,"y":3},{"x":28,"y":2},{"x":29,"y":1},{"x":30,"y":2},{"x":31,"y":2},{"x":32,"y":3},{"x":33,"y":2},{"x":34,"y":2},{"x":35,"y":1},{"x":36,"y":2},{"x":37,"y":2},{"x":38,"y":2},{"x":39,"y":2},{"x":40,"y":2},{"x":41,"y":1},{"x":42,"y":1},{"x":43,"y":1},{"x":44,"y":1},{"x":45,"y":2},{"x":46,"y":1},{"x":47,"y":1},{"x":48,"y":3},{"x":49,"y":2},{"x":50,"y":2},{"x":51,"y":2},{"x":52,"y":2},{"x":53,"y":1},{"x":54,"y":2},{"x":55,"y":1},{"x":56,"y":2},{"x":57,"y":1},{"x":58,"y":1},{"x":59,"y":2},{"x":60,"y":1},{"x":61,"y":2},{"x":62,"y":1},{"x":63,"y":5},{"x":64,"y":3},{"x":65,"y":2},{"x":66,"y":1},{"x":67,"y":1},{"x":68,"y":1},{"x":69,"y":2},{"x":70,"y":1},{"x":71,"y":1},{"x":72,"y":1},{"x":73,"y":2},{"x":74,"y":1},{"x":75,"y":2},{"x":76,"y":1},{"x":77,"y":2},{"x":78,"y":1},{"x":79,"y":2},{"x":80,"y":2},{"x":81,"y":2},{"x":82,"y":1},{"x":83,"y":1},{"x":84,"y":2},{"x":85,"y":2},{"x":86,"y":1},{"x":87,"y":1},{"x":88,"y":2},{"x":89,"y":1},{"x":90,"y":1},{"x":91,"y":1},{"x":92,"y":1},{"x":93,"y":2},{"x":94,"y":1},{"x":95,"y":1},{"x":96,"y":18},{"x":97,"y":2},{"x":98,"y":2},{"x":99,"y":1}]},{"name":"169bae33-b920-41b3-9277-44f1f55a6d71","values":[{"x":0,"y":9850},{"x":1,"y":12514},{"x":2,"y":5306},{"x":3,"y":7484},{"x":4,"y":11165},{"x":5,"y":15874},{"x":6,"y":11526},{"x":7,"y":15863},{"x":8,"y":18114},{"x":9,"y":14407},{"x":10,"y":19217},{"x":11,"y":17995},{"x":12,"y":19799},{"x":13,"y":32958},{"x":14,"y":29486},{"x":15,"y":20248},{"x":16,"y":28116},{"x":17,"y":25913},{"x":18,"y":20442},{"x":19,"y":30220},{"x":20,"y":27254},{"x":21,"y":28558},{"x":22,"y":32759},{"x":23,"y":30340},{"x":24,"y":20987},{"x":25,"y":22044},{"x":26,"y":19185},{"x":27,"y":28354},{"x":28,"y":23745},{"x":29,"y":30834},{"x":30,"y":12223},{"x":31,"y":9205},{"x":32,"y":11352},{"x":33,"y":9800},{"x":34,"y":16039},{"x":35,"y":18602},{"x":36,"y":16760},{"x":37,"y":26537},{"x":38,"y":19986},{"x":39,"y":20402},{"x":40,"y":20962},{"x":41,"y":25970},{"x":42,"y":10903},{"x":43,"y":8603},{"x":44,"y":9974},{"x":45,"y":13454},{"x":46,"y":29657},{"x":47,"y":35514},{"x":48,"y":31135},{"x":49,"y":21519},{"x":50,"y":19781},{"x":51,"y":15225},{"x":52,"y":9065},{"x":53,"y":7779},{"x":54,"y":8316},{"x":55,"y":8608},{"x":56,"y":8811},{"x":57,"y":8609},{"x":58,"y":8458},{"x":59,"y":11023},{"x":60,"y":11625},{"x":61,"y":14437},{"x":62,"y":13134},{"x":63,"y":15832},{"x":64,"y":17405},{"x":65,"y":16983},{"x":66,"y":22651},{"x":67,"y":24817},{"x":68,"y":18802},{"x":69,"y":30581},{"x":70,"y":21962},{"x":71,"y":18011},{"x":72,"y":16313},{"x":73,"y":26544},{"x":74,"y":23302},{"x":75,"y":33422},{"x":76,"y":36819},{"x":77,"y":30542},{"x":78,"y":31692},{"x":79,"y":15508},{"x":80,"y":15779},{"x":81,"y":21486},{"x":82,"y":23010},{"x":83,"y":19949},{"x":84,"y":22164},{"x":85,"y":18030},{"x":86,"y":26872},{"x":87,"y":26458},{"x":88,"y":16734},{"x":89,"y":18334},{"x":90,"y":20674},{"x":91,"y":17600},{"x":92,"y":22764},{"x":93,"y":21183},{"x":94,"y":25166},{"x":95,"y":23549},{"x":96,"y":22995},{"x":97,"y":20815},{"x":98,"y":21966},{"x":99,"y":26594}]}],"marks":[{"type":"line","from":{"data":"50111b8e-ed9e-496f-aced-d45a43d2f4d3"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"yellow"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}},{"type":"line","from":{"data":"66bce705-a4ba-43ed-bc8a-0d9f1aee5652"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"red"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}},{"type":"line","from":{"data":"d15a96e5-b6a5-4601-9760-399ad56d5c56"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"green"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}},{"type":"line","from":{"data":"169bae33-b920-41b3-9277-44f1f55a6d71"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"blue"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}}]},"value":"#gorilla_repl.vega.VegaView{:content {:width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}, :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"50111b8e-ed9e-496f-aced-d45a43d2f4d3\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain [0 36842]}], :axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :data ({:name \"50111b8e-ed9e-496f-aced-d45a43d2f4d3\", :values ({:x 0, :y 10030} {:x 1, :y 12567} {:x 2, :y 5347} {:x 3, :y 7512} {:x 4, :y 11195} {:x 5, :y 15906} {:x 6, :y 11554} {:x 7, :y 15888} {:x 8, :y 18141} {:x 9, :y 14429} {:x 10, :y 19244} {:x 11, :y 18015} {:x 12, :y 19820} {:x 13, :y 32980} {:x 14, :y 29509} {:x 15, :y 20270} {:x 16, :y 28139} {:x 17, :y 25937} {:x 18, :y 20463} {:x 19, :y 30241} {:x 20, :y 27274} {:x 21, :y 28576} {:x 22, :y 32776} {:x 23, :y 30361} {:x 24, :y 21013} {:x 25, :y 22067} {:x 26, :y 19214} {:x 27, :y 28385} {:x 28, :y 23763} {:x 29, :y 30851} {:x 30, :y 12249} {:x 31, :y 9257} {:x 32, :y 11395} {:x 33, :y 9844} {:x 34, :y 16059} {:x 35, :y 18623} {:x 36, :y 16781} {:x 37, :y 26557} {:x 38, :y 20008} {:x 39, :y 20427} {:x 40, :y 21009} {:x 41, :y 26002} {:x 42, :y 10923} {:x 43, :y 8629} {:x 44, :y 10014} {:x 45, :y 13493} {:x 46, :y 29692} {:x 47, :y 35534} {:x 48, :y 31156} {:x 49, :y 21538} {:x 50, :y 19802} {:x 51, :y 15251} {:x 52, :y 9091} {:x 53, :y 7819} {:x 54, :y 8359} {:x 55, :y 8648} {:x 56, :y 8863} {:x 57, :y 8637} {:x 58, :y 8484} {:x 59, :y 11055} {:x 60, :y 11655} {:x 61, :y 14464} {:x 62, :y 13165} {:x 63, :y 15865} {:x 64, :y 17436} {:x 65, :y 17009} {:x 66, :y 22685} {:x 67, :y 24846} {:x 68, :y 18818} {:x 69, :y 30613} {:x 70, :y 21982} {:x 71, :y 18029} {:x 72, :y 16339} {:x 73, :y 26569} {:x 74, :y 23324} {:x 75, :y 33448} {:x 76, :y 36842} {:x 77, :y 30567} {:x 78, :y 31721} {:x 79, :y 15527} {:x 80, :y 15795} {:x 81, :y 21516} {:x 82, :y 23033} {:x 83, :y 19970} {:x 84, :y 22183} {:x 85, :y 18054} {:x 86, :y 26898} {:x 87, :y 26488} {:x 88, :y 16755} {:x 89, :y 18360} {:x 90, :y 20701} {:x 91, :y 17628} {:x 92, :y 22793} {:x 93, :y 21204} {:x 94, :y 25189} {:x 95, :y 23571} {:x 96, :y 23032} {:x 97, :y 20836} {:x 98, :y 21992} {:x 99, :y 26615})} {:name \"66bce705-a4ba-43ed-bc8a-0d9f1aee5652\", :values ({:x 0, :y 64} {:x 1, :y 48} {:x 2, :y 37} {:x 3, :y 26} {:x 4, :y 27} {:x 5, :y 28} {:x 6, :y 26} {:x 7, :y 22} {:x 8, :y 23} {:x 9, :y 20} {:x 10, :y 25} {:x 11, :y 18} {:x 12, :y 19} {:x 13, :y 20} {:x 14, :y 21} {:x 15, :y 20} {:x 16, :y 22} {:x 17, :y 22} {:x 18, :y 18} {:x 19, :y 19} {:x 20, :y 17} {:x 21, :y 16} {:x 22, :y 15} {:x 23, :y 19} {:x 24, :y 24} {:x 25, :y 21} {:x 26, :y 26} {:x 27, :y 28} {:x 28, :y 16} {:x 29, :y 16} {:x 30, :y 24} {:x 31, :y 50} {:x 32, :y 40} {:x 33, :y 42} {:x 34, :y 18} {:x 35, :y 20} {:x 36, :y 19} {:x 37, :y 18} {:x 38, :y 20} {:x 39, :y 23} {:x 40, :y 45} {:x 41, :y 31} {:x 42, :y 19} {:x 43, :y 25} {:x 44, :y 39} {:x 45, :y 37} {:x 46, :y 34} {:x 47, :y 19} {:x 48, :y 18} {:x 49, :y 17} {:x 50, :y 19} {:x 51, :y 24} {:x 52, :y 24} {:x 53, :y 39} {:x 54, :y 41} {:x 55, :y 39} {:x 56, :y 50} {:x 57, :y 27} {:x 58, :y 25} {:x 59, :y 30} {:x 60, :y 29} {:x 61, :y 25} {:x 62, :y 30} {:x 63, :y 28} {:x 64, :y 28} {:x 65, :y 24} {:x 66, :y 33} {:x 67, :y 28} {:x 68, :y 15} {:x 69, :y 30} {:x 70, :y 19} {:x 71, :y 17} {:x 72, :y 25} {:x 73, :y 23} {:x 74, :y 21} {:x 75, :y 24} {:x 76, :y 22} {:x 77, :y 23} {:x 78, :y 28} {:x 79, :y 17} {:x 80, :y 14} {:x 81, :y 28} {:x 82, :y 22} {:x 83, :y 20} {:x 84, :y 17} {:x 85, :y 22} {:x 86, :y 25} {:x 87, :y 29} {:x 88, :y 19} {:x 89, :y 25} {:x 90, :y 26} {:x 91, :y 27} {:x 92, :y 28} {:x 93, :y 19} {:x 94, :y 22} {:x 95, :y 21} {:x 96, :y 19} {:x 97, :y 19} {:x 98, :y 24} {:x 99, :y 20})} {:name \"d15a96e5-b6a5-4601-9760-399ad56d5c56\", :values ({:x 0, :y 116} {:x 1, :y 5} {:x 2, :y 4} {:x 3, :y 2} {:x 4, :y 3} {:x 5, :y 4} {:x 6, :y 2} {:x 7, :y 3} {:x 8, :y 4} {:x 9, :y 2} {:x 10, :y 2} {:x 11, :y 2} {:x 12, :y 2} {:x 13, :y 2} {:x 14, :y 2} {:x 15, :y 2} {:x 16, :y 1} {:x 17, :y 2} {:x 18, :y 3} {:x 19, :y 2} {:x 20, :y 3} {:x 21, :y 2} {:x 22, :y 2} {:x 23, :y 2} {:x 24, :y 2} {:x 25, :y 2} {:x 26, :y 3} {:x 27, :y 3} {:x 28, :y 2} {:x 29, :y 1} {:x 30, :y 2} {:x 31, :y 2} {:x 32, :y 3} {:x 33, :y 2} {:x 34, :y 2} {:x 35, :y 1} {:x 36, :y 2} {:x 37, :y 2} {:x 38, :y 2} {:x 39, :y 2} {:x 40, :y 2} {:x 41, :y 1} {:x 42, :y 1} {:x 43, :y 1} {:x 44, :y 1} {:x 45, :y 2} {:x 46, :y 1} {:x 47, :y 1} {:x 48, :y 3} {:x 49, :y 2} {:x 50, :y 2} {:x 51, :y 2} {:x 52, :y 2} {:x 53, :y 1} {:x 54, :y 2} {:x 55, :y 1} {:x 56, :y 2} {:x 57, :y 1} {:x 58, :y 1} {:x 59, :y 2} {:x 60, :y 1} {:x 61, :y 2} {:x 62, :y 1} {:x 63, :y 5} {:x 64, :y 3} {:x 65, :y 2} {:x 66, :y 1} {:x 67, :y 1} {:x 68, :y 1} {:x 69, :y 2} {:x 70, :y 1} {:x 71, :y 1} {:x 72, :y 1} {:x 73, :y 2} {:x 74, :y 1} {:x 75, :y 2} {:x 76, :y 1} {:x 77, :y 2} {:x 78, :y 1} {:x 79, :y 2} {:x 80, :y 2} {:x 81, :y 2} {:x 82, :y 1} {:x 83, :y 1} {:x 84, :y 2} {:x 85, :y 2} {:x 86, :y 1} {:x 87, :y 1} {:x 88, :y 2} {:x 89, :y 1} {:x 90, :y 1} {:x 91, :y 1} {:x 92, :y 1} {:x 93, :y 2} {:x 94, :y 1} {:x 95, :y 1} {:x 96, :y 18} {:x 97, :y 2} {:x 98, :y 2} {:x 99, :y 1})} {:name \"169bae33-b920-41b3-9277-44f1f55a6d71\", :values ({:x 0, :y 9850} {:x 1, :y 12514} {:x 2, :y 5306} {:x 3, :y 7484} {:x 4, :y 11165} {:x 5, :y 15874} {:x 6, :y 11526} {:x 7, :y 15863} {:x 8, :y 18114} {:x 9, :y 14407} {:x 10, :y 19217} {:x 11, :y 17995} {:x 12, :y 19799} {:x 13, :y 32958} {:x 14, :y 29486} {:x 15, :y 20248} {:x 16, :y 28116} {:x 17, :y 25913} {:x 18, :y 20442} {:x 19, :y 30220} {:x 20, :y 27254} {:x 21, :y 28558} {:x 22, :y 32759} {:x 23, :y 30340} {:x 24, :y 20987} {:x 25, :y 22044} {:x 26, :y 19185} {:x 27, :y 28354} {:x 28, :y 23745} {:x 29, :y 30834} {:x 30, :y 12223} {:x 31, :y 9205} {:x 32, :y 11352} {:x 33, :y 9800} {:x 34, :y 16039} {:x 35, :y 18602} {:x 36, :y 16760} {:x 37, :y 26537} {:x 38, :y 19986} {:x 39, :y 20402} {:x 40, :y 20962} {:x 41, :y 25970} {:x 42, :y 10903} {:x 43, :y 8603} {:x 44, :y 9974} {:x 45, :y 13454} {:x 46, :y 29657} {:x 47, :y 35514} {:x 48, :y 31135} {:x 49, :y 21519} {:x 50, :y 19781} {:x 51, :y 15225} {:x 52, :y 9065} {:x 53, :y 7779} {:x 54, :y 8316} {:x 55, :y 8608} {:x 56, :y 8811} {:x 57, :y 8609} {:x 58, :y 8458} {:x 59, :y 11023} {:x 60, :y 11625} {:x 61, :y 14437} {:x 62, :y 13134} {:x 63, :y 15832} {:x 64, :y 17405} {:x 65, :y 16983} {:x 66, :y 22651} {:x 67, :y 24817} {:x 68, :y 18802} {:x 69, :y 30581} {:x 70, :y 21962} {:x 71, :y 18011} {:x 72, :y 16313} {:x 73, :y 26544} {:x 74, :y 23302} {:x 75, :y 33422} {:x 76, :y 36819} {:x 77, :y 30542} {:x 78, :y 31692} {:x 79, :y 15508} {:x 80, :y 15779} {:x 81, :y 21486} {:x 82, :y 23010} {:x 83, :y 19949} {:x 84, :y 22164} {:x 85, :y 18030} {:x 86, :y 26872} {:x 87, :y 26458} {:x 88, :y 16734} {:x 89, :y 18334} {:x 90, :y 20674} {:x 91, :y 17600} {:x 92, :y 22764} {:x 93, :y 21183} {:x 94, :y 25166} {:x 95, :y 23549} {:x 96, :y 22995} {:x 97, :y 20815} {:x 98, :y 21966} {:x 99, :y 26594})}), :marks ({:type \"line\", :from {:data \"50111b8e-ed9e-496f-aced-d45a43d2f4d3\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"yellow\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}} {:type \"line\", :from {:data \"66bce705-a4ba-43ed-bc8a-0d9f1aee5652\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"red\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}} {:type \"line\", :from {:data \"d15a96e5-b6a5-4601-9760-399ad56d5c56\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"green\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}} {:type \"line\", :from {:data \"169bae33-b920-41b3-9277-44f1f55a6d71\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"blue\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}})}}"}
;; <=

;; **
;;; Spread of complexity in population
;; **

;; @@
(plot/histogram (map :complexity (:rabble result)))
;; @@
;; =>
;;; {"type":"vega","content":{"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"6aa6b167-d1aa-493d-ab08-b3971edba2e4","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":{"data":"6aa6b167-d1aa-493d-ab08-b3971edba2e4","field":"data.y"}}],"marks":[{"type":"line","from":{"data":"6aa6b167-d1aa-493d-ab08-b3971edba2e4"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"interpolate":{"value":"step-before"},"fill":{"value":"steelblue"},"fillOpacity":{"value":0.4},"stroke":{"value":"steelblue"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}}],"data":[{"name":"6aa6b167-d1aa-493d-ab08-b3971edba2e4","values":[{"x":1.0,"y":0},{"x":3.4285714285714293,"y":8.0},{"x":5.8571428571428585,"y":16.0},{"x":8.285714285714288,"y":6.0},{"x":10.714285714285717,"y":8.0},{"x":13.142857142857146,"y":4.0},{"x":15.571428571428575,"y":3.0},{"x":18.000000000000004,"y":5.0},{"x":20.428571428571434,"y":0}]}],"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50}},"value":"#gorilla_repl.vega.VegaView{:content {:axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"6aa6b167-d1aa-493d-ab08-b3971edba2e4\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain {:data \"6aa6b167-d1aa-493d-ab08-b3971edba2e4\", :field \"data.y\"}}], :marks [{:type \"line\", :from {:data \"6aa6b167-d1aa-493d-ab08-b3971edba2e4\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :interpolate {:value \"step-before\"}, :fill {:value \"steelblue\"}, :fillOpacity {:value 0.4}, :stroke {:value \"steelblue\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}}], :data [{:name \"6aa6b167-d1aa-493d-ab08-b3971edba2e4\", :values ({:x 1.0, :y 0} {:x 3.4285714285714293, :y 8.0} {:x 5.8571428571428585, :y 16.0} {:x 8.285714285714288, :y 6.0} {:x 10.714285714285717, :y 8.0} {:x 13.142857142857146, :y 4.0} {:x 15.571428571428575, :y 3.0} {:x 18.000000000000004, :y 5.0} {:x 20.428571428571434, :y 0})}], :width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}}}"}
;; <=

;; @@
(plot/list-plot (:mean (:complexity @metrics/metrics)) :joined true)
;; @@
;; =>
;;; {"type":"vega","content":{"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"79879f9c-758a-49f4-884e-b18b78961009","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":{"data":"79879f9c-758a-49f4-884e-b18b78961009","field":"data.y"}}],"marks":[{"type":"line","from":{"data":"79879f9c-758a-49f4-884e-b18b78961009"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"#FF29D2"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}}],"data":[{"name":"79879f9c-758a-49f4-884e-b18b78961009","values":[{"x":0,"y":4.29},{"x":1,"y":4.86},{"x":2,"y":2.71},{"x":3,"y":3.53},{"x":4,"y":4.38},{"x":5,"y":5.71},{"x":6,"y":4.69},{"x":7,"y":5.76},{"x":8,"y":5.54},{"x":9,"y":5.15},{"x":10,"y":5.89},{"x":11,"y":5.53},{"x":12,"y":5.68},{"x":13,"y":8.67},{"x":14,"y":7.66},{"x":15,"y":6.36},{"x":16,"y":7.66},{"x":17,"y":7.37},{"x":18,"y":6.11},{"x":19,"y":7.86},{"x":20,"y":6.68},{"x":21,"y":7.67},{"x":22,"y":8.01},{"x":23,"y":7.75},{"x":24,"y":6.49},{"x":25,"y":6.0},{"x":26,"y":5.67},{"x":27,"y":6.9},{"x":28,"y":5.7},{"x":29,"y":6.6},{"x":30,"y":3.93},{"x":31,"y":3.31},{"x":32,"y":3.81},{"x":33,"y":3.72},{"x":34,"y":5.6},{"x":35,"y":5.59},{"x":36,"y":5.12},{"x":37,"y":6.32},{"x":38,"y":5.58},{"x":39,"y":5.05},{"x":40,"y":5.61},{"x":41,"y":6.08},{"x":42,"y":3.88},{"x":43,"y":3.47},{"x":44,"y":3.78},{"x":45,"y":4.48},{"x":46,"y":6.35},{"x":47,"y":7.45},{"x":48,"y":6.45},{"x":49,"y":4.88},{"x":50,"y":5.14},{"x":51,"y":4.5},{"x":52,"y":3.81},{"x":53,"y":3.44},{"x":54,"y":3.6},{"x":55,"y":3.76},{"x":56,"y":3.67},{"x":57,"y":3.61},{"x":58,"y":3.75},{"x":59,"y":4.38},{"x":60,"y":4.19},{"x":61,"y":4.79},{"x":62,"y":4.67},{"x":63,"y":4.85},{"x":64,"y":5.38},{"x":65,"y":5.04},{"x":66,"y":5.27},{"x":67,"y":5.85},{"x":68,"y":4.91},{"x":69,"y":6.57},{"x":70,"y":5.46},{"x":71,"y":5.19},{"x":72,"y":4.77},{"x":73,"y":6.28},{"x":74,"y":5.38},{"x":75,"y":6.44},{"x":76,"y":6.97},{"x":77,"y":6.13},{"x":78,"y":6.55},{"x":79,"y":5.0},{"x":80,"y":4.91},{"x":81,"y":5.85},{"x":82,"y":5.99},{"x":83,"y":5.45},{"x":84,"y":5.86},{"x":85,"y":5.47},{"x":86,"y":6.35},{"x":87,"y":6.25},{"x":88,"y":5.47},{"x":89,"y":5.72},{"x":90,"y":6.24},{"x":91,"y":5.78},{"x":92,"y":6.5},{"x":93,"y":6.06},{"x":94,"y":6.79},{"x":95,"y":6.32},{"x":96,"y":6.35},{"x":97,"y":5.78},{"x":98,"y":6.27},{"x":99,"y":6.55}]}],"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50}},"value":"#gorilla_repl.vega.VegaView{:content {:axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"79879f9c-758a-49f4-884e-b18b78961009\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain {:data \"79879f9c-758a-49f4-884e-b18b78961009\", :field \"data.y\"}}], :marks [{:type \"line\", :from {:data \"79879f9c-758a-49f4-884e-b18b78961009\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"#FF29D2\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}}], :data [{:name \"79879f9c-758a-49f4-884e-b18b78961009\", :values ({:x 0, :y 4.29} {:x 1, :y 4.86} {:x 2, :y 2.71} {:x 3, :y 3.53} {:x 4, :y 4.38} {:x 5, :y 5.71} {:x 6, :y 4.69} {:x 7, :y 5.76} {:x 8, :y 5.54} {:x 9, :y 5.15} {:x 10, :y 5.89} {:x 11, :y 5.53} {:x 12, :y 5.68} {:x 13, :y 8.67} {:x 14, :y 7.66} {:x 15, :y 6.36} {:x 16, :y 7.66} {:x 17, :y 7.37} {:x 18, :y 6.11} {:x 19, :y 7.86} {:x 20, :y 6.68} {:x 21, :y 7.67} {:x 22, :y 8.01} {:x 23, :y 7.75} {:x 24, :y 6.49} {:x 25, :y 6.0} {:x 26, :y 5.67} {:x 27, :y 6.9} {:x 28, :y 5.7} {:x 29, :y 6.6} {:x 30, :y 3.93} {:x 31, :y 3.31} {:x 32, :y 3.81} {:x 33, :y 3.72} {:x 34, :y 5.6} {:x 35, :y 5.59} {:x 36, :y 5.12} {:x 37, :y 6.32} {:x 38, :y 5.58} {:x 39, :y 5.05} {:x 40, :y 5.61} {:x 41, :y 6.08} {:x 42, :y 3.88} {:x 43, :y 3.47} {:x 44, :y 3.78} {:x 45, :y 4.48} {:x 46, :y 6.35} {:x 47, :y 7.45} {:x 48, :y 6.45} {:x 49, :y 4.88} {:x 50, :y 5.14} {:x 51, :y 4.5} {:x 52, :y 3.81} {:x 53, :y 3.44} {:x 54, :y 3.6} {:x 55, :y 3.76} {:x 56, :y 3.67} {:x 57, :y 3.61} {:x 58, :y 3.75} {:x 59, :y 4.38} {:x 60, :y 4.19} {:x 61, :y 4.79} {:x 62, :y 4.67} {:x 63, :y 4.85} {:x 64, :y 5.38} {:x 65, :y 5.04} {:x 66, :y 5.27} {:x 67, :y 5.85} {:x 68, :y 4.91} {:x 69, :y 6.57} {:x 70, :y 5.46} {:x 71, :y 5.19} {:x 72, :y 4.77} {:x 73, :y 6.28} {:x 74, :y 5.38} {:x 75, :y 6.44} {:x 76, :y 6.97} {:x 77, :y 6.13} {:x 78, :y 6.55} {:x 79, :y 5.0} {:x 80, :y 4.91} {:x 81, :y 5.85} {:x 82, :y 5.99} {:x 83, :y 5.45} {:x 84, :y 5.86} {:x 85, :y 5.47} {:x 86, :y 6.35} {:x 87, :y 6.25} {:x 88, :y 5.47} {:x 89, :y 5.72} {:x 90, :y 6.24} {:x 91, :y 5.78} {:x 92, :y 6.5} {:x 93, :y 6.06} {:x 94, :y 6.79} {:x 95, :y 6.32} {:x 96, :y 6.35} {:x 97, :y 5.78} {:x 98, :y 6.27} {:x 99, :y 6.55})}], :width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}}}"}
;; <=

;; @@

;; @@
