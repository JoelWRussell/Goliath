;; gorilla-repl.fileformat = 1

;; **
;;; #Networked Hybrid Windows 2
;;; 
;;; 
;; **

;; @@
(def experimentalDataSz "resources/mma_double.csv")
(def dt 0.1)
(def df 2)

;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/df</span>","value":"#'user/df"}
;; <=

;; @@
(import '[com.lagrangianmining.GorillaClient])
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@
(def client (com.lagrangianmining.GorillaClient.))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/client</span>","value":"#'user/client"}
;; <=

;; @@
(defn networkScore [gen]
(let[
sizePoly (map #(count (flatten %)) gen)]
(.NewZeitgeist user/client (count gen) (into-array Integer/TYPE sizePoly) (into-array Integer/TYPE (flatten gen)) 2 0.1)
(into [] (.GetScores user/client))
)
)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/networkScore</span>","value":"#'user/networkScore"}
;; <=

;; @@
(.InitMathKernel client)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>true</span>","value":"true"}
;; <=

;; @@
(.NewData client "resources/mma_double.csv")
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>true</span>","value":"true"}
;; <=

;; @@
(.PrepareData client false df dt)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>true</span>","value":"true"}
;; <=

;; @@
(require '[darwin.evolution.core :as evolution]
         '[darwin.evolution.metrics :as metrics]
         '[darwin.evolution.reproduction :as reproduction]
         '[darwin.evolution.scoring :as scoring]
         '[darwin.evolution.selection :as selection]
         '[darwin.evolution.transform :as transform]
         '[darwin.evolution.pareto :as pareto]
         '[darwin.algorithms.spea2 :as spea2]
    	 '[darwin.core :as darwin]
         '[darwin.evolution.metrics :as metrics]
         '[algebolic.expression.core :as expression]
         '[algebolic.expression.tree :as tree]
         '[algebolic.expression.genetics :as genetics]
         '[algebolic.expression.score :as score]
         '[algebolic.expression.render :as spea-render]
         '[algebolic.expression.interpreter :as interpreter]
         '[criterium.core :as criterium]
         )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@
(def prob_inheritance 0.75) ;;to do with crossover  ;probability that terms are directly inherited from parent to child.
(def mutate_pref1 0.8)
(def mutate_pref2 0.8);;change 1 df within a gene  ;prob. that "gene-replace"/"gene-add" is called when mutate is called.
(def mutate_pref3 0.2) ;; add a whole new gene/replace a term (big change)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/mutate_pref3</span>","value":"#'user/mutate_pref3"}
;; <=

;; @@
(def power_sum 10)    ;the maximum sum of powers present in each term within a polynomial.
(def poly_length 10)  ;the maximum number of terms of a polynomial. 
(def power 8) 		  ;the maximum power a variable can be raised to.
(def term_length (* 2 user/df))   ;the number of variables within each term - 2*degree_freedom
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/term_length</span>","value":"#'user/term_length"}
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
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/new-poly-term</span>","value":"#'user/new-poly-term"}
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
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/create-genotype</span>","value":"#'user/create-genotype"}
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
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/distribute-parent</span>","value":"#'user/distribute-parent"}
;; <=

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
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/cross-over</span>","value":"#'user/cross-over"}
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
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/gene-replace</span>","value":"#'user/gene-replace"}
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
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/gene-add</span>","value":"#'user/gene-add"}
;; <=

;; @@
(defn gene-tweak
  [indv]
  (let [n (rand-int (count indv)) gene (nth indv n) new-gene (assoc gene (rand-int (count gene)) (rand-int power))]
    (assoc indv n new-gene)
    ))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/gene-tweak</span>","value":"#'user/gene-tweak"}
;; <=

;; @@
(defn gene-delete [indv] 
  (let [new-indv (subvec indv 1 (count indv))]
    (if (= new-indv []) indv
      new-indv))
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/gene-delete</span>","value":"#'user/gene-delete"}
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
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/mutate</span>","value":"#'user/mutate"}
;; <=

;; @@
(defn random-initial-population 
  
   "creates an innitial population of polynimials
  size - number of individuals in the population"
   
  [size max_length max_power max_sum_powers num_terms]
 
  (repeatedly size #(create-genotype max_length max_power max_sum_powers num_terms)) 
       
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/random-initial-population</span>","value":"#'user/random-initial-population"}
;; <=

;; @@
(def initial-zeitgeist (evolution/make-zeitgeist (random-initial-population 100 poly_length power power_sum term_length)))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/initial-zeitgeist</span>","value":"#'user/initial-zeitgeist"}
;; <=

;; @@
(def generation-config			      
  (let [size-limit 120                   
        min-size 20
        ea-config (spea2/spea2-config
                    {:goals [:error :complexity]
                     :archive-size 100
                     :binary-ops [{:op cross-over :repeat 40}]
                     :unary-ops [{:op mutate :repeat 20}]})
        score-functions {:complexity (fn [x] (count x))
                         :error (fn [e] (* 1.0 e))}]
    {:ea-config              ea-config
     :score-functions        score-functions
     :reporting-function     (fn [z] (print ".") (flush))}))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/generation-config</span>","value":"#'user/generation-config"}
;; <=

;; @@
(time (def result (evolution/run-evolution generation-config initial-zeitgeist (fn [zg gc] (>= (:age zg) 1)))))
;; @@
;; ->
;;; ({:error 0.028232322722737884, :complexity 12, :genotype [[0 0 0 1] [0 0 1 0] [0 1 0 0]]} {:error -0.04899291061013166, :complexity 20, :genotype [[0 0 1 0] [1 0 0 0] [0 0 0 1] [1 0 6 0] [2 5 2 0]]} {:error -1.4399696055586102, :complexity 20, :genotype [[0 0 0 1] [0 0 1 0] [0 1 0 0] [1 0 0 0] [3 1 0 0]]} {:error -2.925055521072637E-4, :complexity 12, :genotype [[0 2 7 1] [0 0 1 0] [1 0 0 0]]} {:error 0.0, :complexity 8, :genotype [[1 0 0 0] [0 1 0 0]]} {:error 28.211219304298957, :complexity 12, :genotype [[0 0 1 0] [0 0 0 1] [4 1 3 2]]} {:error -0.016436085323039997, :complexity 28, :genotype [[1 4 0 1] [0 0 1 0] [5 0 3 0] [0 0 4 0] [1 0 0 0] [0 0 0 1] [0 1 0 0]]} {:error -7.733736949462403E-4, :complexity 16, :genotype [[3 0 0 7] [1 0 0 0] [0 0 1 0] [0 1 0 0]]} {:error -6.163854529145253E-4, :complexity 20, :genotype [[0 0 0 1] [0 5 1 1] [1 0 0 0] [2 0 6 0] [2 0 4 2]]} {:error -0.01310161150949659, :complexity 24, :genotype [[0 0 0 1] [1 0 0 0] [0 1 0 0] [0 5 4 1] [0 0 1 0] [1 2 0 6]]} {:error -0.003744219084468871, :complexity 20, :genotype [[2 4 0 1] [0 0 1 0] [1 0 0 0] [0 1 0 0] [0 0 0 1]]} {:error -0.004218235824945162, :complexity 16, :genotype [[1 0 3 0] [0 1 0 0] [1 0 0 0] [1 3 0 2]]} {:error -0.023512378913830134, :complexity 28, :genotype [[2 3 3 2] [2 3 0 5] [0 0 0 1] [0 1 7 0] [4 3 1 1] [0 1 0 0] [2 6 0 2]]} {:error 0.0, :complexity 8, :genotype [[1 0 0 0] [0 1 0 0]]} {:error 0.0, :complexity 8, :genotype [[1 0 0 0] [0 1 0 0]]} {:error 0.028232322722737884, :complexity 12, :genotype [[0 0 1 0] [0 1 0 0] [0 0 0 1]]} {:error 18.502165076507065, :complexity 8, :genotype [[0 0 1 0] [0 5 3 0]]} {:error 16.71784908976807, :complexity 8, :genotype [[0 0 0 1] [7 0 2 1]]} {:error 1.5832220019188459, :complexity 4, :genotype [[1 0 0 0]]} {:error -0.0025642878594982982, :complexity 20, :genotype [[0 0 0 1] [0 2 6 0] [0 0 1 0] [0 1 0 0] [1 1 4 3]]} {:error 0.00986814613767775, :complexity 8, :genotype [[0 0 0 1] [0 1 0 0]]} {:error -1.8146670833770209, :complexity 24, :genotype [[0 0 0 1] [0 0 1 0] [2 5 2 1] [0 1 0 0] [2 3 0 0] [0 3 1 5]]} {:error 0.028232322722737884, :complexity 12, :genotype [[0 0 1 0] [0 1 0 0] [0 0 0 1]]} {:error -0.03735029749918881, :complexity 16, :genotype [[0 1 0 0] [1 0 0 0] [2 3 0 4] [0 0 1 0]]} {:error 0.00986814613767775, :complexity 8, :genotype [[1 0 0 0] [0 0 1 0]]} {:error -2.9690222798052447E-4, :complexity 20, :genotype [[0 1 0 0] [1 0 3 1] [0 6 1 0] [3 3 3 1] [0 0 0 1]]} {:error 0.0, :complexity 16, :genotype [[0 0 1 0] [0 1 0 0] [1 0 0 0] [0 0 0 1]]} {:error -0.004518896012615737, :complexity 20, :genotype [[1 0 0 0] [0 0 1 0] [0 1 0 0] [0 0 0 1] [0 5 3 0]]} {:error 1.5832220019188459, :complexity 4, :genotype [[1 0 0 0]]} {:error -0.001369312688205371, :complexity 16, :genotype [[0 0 0 1] [1 0 0 0] [0 2 6 1] [0 0 1 0]]} {:error 0.0, :complexity 8, :genotype [[0 1 0 0] [1 0 0 0]]} {:error -6.355407162029653E-4, :complexity 12, :genotype [[1 0 0 0] [1 0 7 2] [0 1 0 0]]} {:error 0.0, :complexity 12, :genotype [[0 0 1 0] [1 0 0 0] [0 1 0 0]]} {:error -0.013615412403065573, :complexity 20, :genotype [[0 0 0 1] [0 1 0 0] [1 0 0 0] [0 0 4 3] [0 0 1 0]]} {:error -0.010280170393330347, :complexity 20, :genotype [[0 0 0 1] [0 0 1 0] [1 0 5 4] [1 0 0 0] [4 0 0 6]]} {:error 16.212236546842988, :complexity 4, :genotype [[0 0 1 0]]} {:error 0.0, :complexity 16, :genotype [[0 0 0 1] [0 0 1 0] [1 0 0 0] [0 1 0 0]]} {:error 0.0, :complexity 16, :genotype [[1 0 0 0] [0 1 0 0] [0 0 0 1] [0 0 1 0]]} {:error -0.10715149932514655, :complexity 16, :genotype [[1 0 0 0] [0 0 1 0] [1 6 0 2] [0 1 0 0]]} {:error -4.093642274742216E-4, :complexity 28, :genotype [[0 1 0 0] [3 3 3 1] [0 0 2 2] [0 0 1 0] [1 1 6 2] [1 0 0 0] [0 0 0 1]]} {:error 0.0, :complexity 16, :genotype [[1 0 0 0] [0 0 1 0] [0 0 0 1] [0 1 0 0]]} {:error 0.00986814613767775, :complexity 8, :genotype [[1 0 0 0] [0 0 1 0]]} {:error 0.0, :complexity 12, :genotype [[1 0 0 0] [0 0 1 0] [0 1 0 0]]} {:error 1.5832220019188459, :complexity 4, :genotype [[1 0 0 0]]} {:error 16.212236546842988, :complexity 4, :genotype [[0 0 1 0]]} {:error -0.0016329927732507256, :complexity 12, :genotype [[1 0 0 0] [2 0 0 1] [4 2 0 2]]} {:error -0.03014074658160298, :complexity 16, :genotype [[0 1 0 0] [0 0 0 1] [0 3 4 0] [0 0 1 0]]} {:error -0.0011123062761258483, :complexity 20, :genotype [[0 0 1 0] [0 1 0 0] [1 1 4 2] [1 0 0 0] [0 0 0 1]]} {:error 1.5832220019188459, :complexity 4, :genotype [[1 0 0 0]]} {:error -0.006806603547829176, :complexity 20, :genotype [[0 0 1 0] [7 0 0 3] [1 0 0 0] [0 0 0 1] [0 1 0 0]]} {:error 11.565456733340245, :complexity 12, :genotype [[0 0 1 0] [0 4 2 4] [1 5 0 2]]} {:error -0.0017870704353952817, :complexity 20, :genotype [[0 1 0 0] [0 0 0 1] [0 0 1 0] [1 0 0 0] [4 3 3 0]]} {:error 0.0, :complexity 8, :genotype [[0 1 0 0] [1 0 0 0]]} {:error 0.028232322722737884, :complexity 12, :genotype [[0 0 1 0] [0 0 0 1] [0 1 0 0]]} {:error 0.0, :complexity 16, :genotype [[0 1 0 0] [0 0 0 1] [0 0 1 0] [1 0 0 0]]} {:error 0.0, :complexity 16, :genotype [[1 0 0 0] [0 0 0 1] [0 0 1 0] [0 1 0 0]]} {:error -0.004583687484335895, :complexity 12, :genotype [[2 3 0 3] [1 0 0 0] [0 1 0 0]]} {:error -0.8765963074677557, :complexity 24, :genotype [[0 0 0 1] [0 0 1 0] [4 0 1 1] [1 0 0 0] [2 4 2 1] [6 0 0 0]]} {:error 0.00986814613767775, :complexity 8, :genotype [[1 0 0 0] [0 0 0 1]]} {:error -3.4282408092324927, :complexity 20, :genotype [[1 0 0 0] [0 1 0 0] [0 3 0 4] [0 0 1 0] [3 5 0 0]]} {:error -0.0017427417335262538, :complexity 16, :genotype [[5 1 1 1] [0 0 1 0] [0 1 2 5] [1 0 0 0]]} {:error -4.355690086403119E-4, :complexity 20, :genotype [[1 0 0 0] [0 4 2 4] [1 5 1 3] [0 0 1 0] [0 2 2 4]]} {:error 0.028232322722737884, :complexity 12, :genotype [[0 0 1 0] [0 0 0 1] [1 0 0 0]]} {:error -0.0014438502994683718, :complexity 20, :genotype [[0 0 1 0] [0 0 0 1] [4 0 2 4] [1 0 0 0] [6 0 4 0]]} {:error -1.7613582620178936, :complexity 20, :genotype [[1 3 3 3] [6 1 0 0] [0 0 1 0] [0 1 0 0] [0 0 0 1]]} {:error -0.0747173844611858, :complexity 16, :genotype [[0 0 1 0] [0 0 0 1] [1 6 0 2] [0 1 0 0]]} {:error 1.5832220019188459, :complexity 4, :genotype [[0 1 0 0]]} {:error -9.722469016059903E-6, :complexity 16, :genotype [[0 0 0 1] [0 1 0 0] [0 0 1 0] [1 5 2 0]]} {:error 37.771673155871085, :complexity 16, :genotype [[1 1 4 2] [1 0 8 0] [0 0 0 1] [0 0 1 0]]} {:error 0.028232322722737884, :complexity 12, :genotype [[1 0 0 0] [0 0 1 0] [0 0 0 1]]} {:error -1.2848993246911792E-4, :complexity 12, :genotype [[0 1 0 0] [1 0 0 0] [3 1 4 2]]} {:error 0.0, :complexity 12, :genotype [[0 1 0 0] [3 0 4 0] [0 0 1 0]]} {:error 16.212236546842988, :complexity 4, :genotype [[0 0 1 0]]} {:error -0.012524971127530439, :complexity 16, :genotype [[0 1 0 0] [1 0 0 0] [1 0 0 6] [0 0 0 1]]} {:error -5.217282730569098E-4, :complexity 16, :genotype [[0 4 1 2] [1 0 0 0] [0 0 0 1] [0 1 0 0]]} {:error 0.0, :complexity 12, :genotype [[0 1 0 0] [0 0 1 0] [1 0 0 0]]} {:error 0.0, :complexity 16, :genotype [[1 0 0 0] [0 0 1 0] [0 0 0 1] [0 1 0 0]]} {:error 16.212236546842988, :complexity 4, :genotype [[0 0 1 0]]} {:error 0.0, :complexity 16, :genotype [[0 0 0 1] [1 0 0 0] [0 1 0 0] [0 0 1 0]]} {:error 0.0, :complexity 12, :genotype [[1 0 0 0] [0 0 0 1] [0 1 0 0]]} {:error 0.0, :complexity 8, :genotype [[1 0 0 0] [0 1 0 0]]} {:error 33.81130242328212, :complexity 12, :genotype [[0 0 0 1] [0 0 1 0] [2 0 5 0]]} {:error 0.0, :complexity 16, :genotype [[0 0 1 0] [0 0 0 1] [0 1 0 0] [1 0 0 0]]} {:error -0.004564171439618814, :complexity 20, :genotype [[6 0 2 2] [0 1 0 0] [1 0 0 0] [0 0 1 0] [0 0 0 1]]} {:error 0.0, :complexity 12, :genotype [[0 1 0 0] [1 0 0 0] [0 0 1 0]]} {:error -0.010450184383633588, :complexity 12, :genotype [[0 3 6 0] [1 0 0 0] [0 1 0 0]]} {:error 0.028232322722737884, :complexity 12, :genotype [[0 0 1 0] [0 1 0 0] [0 0 0 1]]} {:error 16.212236546842988, :complexity 4, :genotype [[0 0 0 1]]} {:error 0.028232322722737884, :complexity 12, :genotype [[0 0 1 0] [0 0 0 1] [0 1 0 0]]} {:error 16.212236546842988, :complexity 4, :genotype [[0 0 1 0]]} {:error -0.005988672929277489, :complexity 20, :genotype [[0 0 0 1] [1 0 0 0] [6 0 0 4] [2 0 0 7] [0 1 0 0]]} {:error 0.0, :complexity 16, :genotype [[1 0 0 0] [0 1 0 0] [0 0 1 0] [0 0 0 1]]} {:error -0.004129929703584058, :complexity 16, :genotype [[0 0 1 0] [1 0 0 0] [1 1 3 4] [5 2 1 0]]} {:error -4.614989581352499E-5, :complexity 16, :genotype [[1 0 0 0] [0 1 0 0] [1 3 3 3] [0 0 0 1]]} {:error 4.005377030182627, :complexity 12, :genotype [[1 5 2 0] [3 5 0 1] [0 0 0 1]]} {:error -4.1096567855760434E-5, :complexity 16, :genotype [[0 2 1 7] [0 0 1 0] [0 1 0 0] [0 0 0 1]]} {:error 16.212236546842988, :complexity 4, :genotype [[0 0 1 0]]} {:error 1.5832220019188459, :complexity 4, :genotype [[1 0 0 0]]} {:error 0.0, :complexity 12, :genotype [[0 0 1 0] [0 1 0 0] [1 0 0 0]]} {:error 16.212236546842988, :complexity 8, :genotype [[0 0 1 0] [0 0 0 1]]})
;;; ({:error 0.00986814613767775, :complexity 8, :genotype [[0 0 0 1] [1 0 0 0]], :age 1} {:error -0.030141422585557676, :complexity 16, :genotype [[0 1 0 0] [1 0 0 0] [0 3 4 0] [0 0 1 0]], :age 1} {:error -1.4399696055586102, :complexity 20, :genotype [[0 0 0 1] [0 0 1 0] [0 1 0 0] [1 0 0 0] [3 1 0 0]], :age 1} {:error -1.7087900182704232E-4, :complexity 12, :genotype [[0 1 0 0] [0 0 0 1] [4 1 3 2]], :age 1} {:error -4.813730779767095E-5, :complexity 12, :genotype [[0 0 1 0] [0 1 2 5] [1 0 0 0]], :age 1} {:error 1.5832220019188459, :complexity 4, :genotype [[1 0 0 0]], :age 1} {:error 0.028232322722737884, :complexity 12, :genotype [[0 1 0 0] [0 0 1 0] [0 0 0 1]], :age 1} {:error -4.813730779767095E-5, :complexity 16, :genotype [[0 0 0 1] [0 0 1 0] [0 1 2 5] [1 0 0 0]], :age 1} {:error -0.0016931323696107634, :complexity 16, :genotype [[1 0 3 0] [5 2 1 2] [1 0 0 0] [1 3 0 2]], :age 1} {:error 1.5832220019188459, :complexity 4, :genotype [[1 0 0 0]], :age 1} {:error 0.0, :complexity 12, :genotype [[1 0 0 0] [0 0 0 1] [0 1 0 0]], :age 1} {:error 0.028232322722737884, :complexity 12, :genotype [[0 0 0 1] [0 0 1 0] [1 0 0 0]], :age 1} {:error -0.01310161150949659, :complexity 24, :genotype [[0 0 0 1] [1 0 0 0] [0 1 0 0] [0 5 4 1] [0 0 1 0] [1 2 0 6]], :age 1} {:error -0.0463274789387608, :complexity 16, :genotype [[0 1 0 0] [3 2 0 1] [2 3 0 4] [0 0 1 0]], :age 1} {:error 1.5832220019188459, :complexity 4, :genotype [[1 0 0 0]], :age 1} {:error 0.0, :complexity 16, :genotype [[0 0 1 0] [1 0 0 0] [0 1 0 0] [0 0 0 1]], :age 1} {:error 16.212236546842988, :complexity 12, :genotype [[3 0 5 2] [0 0 1 0] [0 0 0 1]], :age 1} {:error -0.02349870310965491, :complexity 24, :genotype [[2 3 0 5] [0 0 0 1] [0 1 7 0] [4 3 1 1] [0 1 0 0] [2 6 0 2]], :age 1} {:error -0.001609254832621151, :complexity 16, :genotype [[1 0 0 0] [0 0 1 0] [2 2 2 4] [0 1 0 0]], :age 1} {:error 0.0, :complexity 12, :genotype [[1 0 0 0] [0 1 0 0] [0 0 1 0]], :age 1} {:error 0.00986814613767775, :complexity 8, :genotype [[1 0 0 0] [0 0 1 0]], :age 1} {:error 1.5832220019188459, :complexity 4, :genotype [[0 1 0 0]], :age 1} {:error -1.1788791374362439E-4, :complexity 8, :genotype [[1 0 0 0] [3 1 4 2]], :age 1} {:error -0.0025642878594982982, :complexity 20, :genotype [[0 1 0 0] [0 0 1 0] [0 2 6 0] [0 0 0 1] [1 1 4 3]], :age 1} {:error 0.0, :complexity 8, :genotype [[1 0 0 0] [0 1 0 0]], :age 1} {:error 0.0, :complexity 8, :genotype [[1 0 0 0] [0 1 0 0]], :age 1} {:error 0.0, :complexity 8, :genotype [[1 0 0 0] [0 1 0 0]], :age 1} {:error 0.0, :complexity 8, :genotype [[1 0 0 0] [0 1 0 0]], :age 1} {:error 16.212236546842988, :complexity 4, :genotype [[0 0 1 0]], :age 1} {:error 0.0, :complexity 8, :genotype [[0 1 0 0] [1 0 0 0]], :age 1} {:error 0.00986814613767775, :complexity 8, :genotype [[1 0 0 0] [0 0 1 0]], :age 1} {:error 16.212236546842988, :complexity 8, :genotype [[0 2 6 1] [0 0 0 1]], :age 1} {:error -0.0037025915022585914, :complexity 12, :genotype [[0 1 0 0] [2 3 0 3] [0 0 0 1]], :age 1} {:error 1.5832220019188459, :complexity 4, :genotype [[1 0 0 0]], :age 1} {:error 1.5832220019188459, :complexity 4, :genotype [[0 1 0 0]], :age 1} {:error 0.0, :complexity 8, :genotype [[1 0 0 0] [0 1 0 0]], :age 1} {:error -4.614989581352499E-5, :complexity 16, :genotype [[1 0 0 0] [0 1 0 0] [1 3 3 3] [0 0 0 1]], :age 1} {:error -0.001443850299468483, :complexity 16, :genotype [[1 0 0 0] [0 0 1 0] [6 0 4 0] [4 0 2 4]], :age 1} {:error 16.212236546842988, :complexity 4, :genotype [[0 0 1 0]], :age 1} {:error 1.5832220019188459, :complexity 4, :genotype [[1 0 0 0]], :age 1} {:error -0.10715149932514655, :complexity 16, :genotype [[1 0 0 0] [0 1 0 0] [1 6 0 2] [0 0 1 0]], :age 1} {:error 0.0, :complexity 12, :genotype [[1 0 0 0] [0 1 0 0] [0 0 1 0]], :age 1} {:error -0.010280170393330235, :complexity 16, :genotype [[1 0 0 0] [0 0 1 0] [4 0 0 6] [1 0 5 4]], :age 1} {:error 0.0, :complexity 12, :genotype [[1 0 0 0] [0 1 0 0] [0 0 0 1]], :age 1} {:error 0.0, :complexity 8, :genotype [[1 0 0 0] [0 1 0 0]], :age 1} {:error -1.1788791374362439E-4, :complexity 8, :genotype [[1 0 0 0] [3 1 4 2]], :age 1} {:error 0.0, :complexity 12, :genotype [[1 0 0 0] [0 1 0 0] [0 0 1 0]], :age 1} {:error 25.877569097487626, :complexity 20, :genotype [[0 0 1 0] [4 0 0 6] [1 0 5 4] [0 0 0 1] [0 3 4 0]], :age 1} {:error -0.0010352610931476571, :complexity 12, :genotype [[0 1 0 0] [0 1 7 0] [2 3 0 5]], :age 1} {:error -0.05832939462142964, :complexity 24, :genotype [[1 0 0 0] [0 0 1 0] [2 3 3 2] [0 0 0 1] [4 3 1 1] [2 6 0 2]], :age 1} {:error 0.0, :complexity 12, :genotype [[1 0 0 0] [0 1 0 0] [0 0 1 0]], :age 1} {:error 11.655285311541613, :complexity 4, :genotype [[1 6 0 2]], :age 1} {:error 12.23215435926609, :complexity 8, :genotype [[0 0 1 0] [4 0 1 1]], :age 1} {:error -0.7866147014246379, :complexity 24, :genotype [[2 4 2 1] [1 0 0 0] [0 1 0 0] [0 0 1 0] [6 0 0 0] [0 0 0 1]], :age 1} {:error -1.7636524410673601E-4, :complexity 8, :genotype [[1 0 0 0] [5 2 1 0]], :age 1} {:error -0.007746184793130387, :complexity 16, :genotype [[1 0 0 0] [0 1 0 0] [0 0 1 0] [1 1 3 4]], :age 1} {:error 0.0, :complexity 12, :genotype [[1 0 0 0] [0 1 0 0] [0 0 1 0]], :age 1} {:error 1.5832220019188459, :complexity 4, :genotype [[1 0 0 0]], :age 1} {:error 0.0, :complexity 8, :genotype [[1 0 0 0] [0 1 0 0]], :age 1} {:error 0.0, :complexity 8, :genotype [[1 0 0 0] [0 1 0 0]], :age 1} {:error 0.00986814613767775, :complexity 8, :genotype [[1 0 0 0] [0 0 1 0]], :age 1} {:error 0.0, :complexity 8, :genotype [[1 0 0 0] [0 1 0 0]], :age 1} {:error 0.0, :complexity 8, :genotype [[1 0 0 0] [0 1 0 0]], :age 1} {:error -0.005487134652710372, :complexity 20, :genotype [[0 1 0 0] [6 0 0 4] [2 0 0 7] [1 0 7 2] [0 0 0 1]], :age 1} {:error 9.360885083617099, :complexity 12, :genotype [[0 1 7 0] [2 3 3 2] [2 6 0 2]], :age 1} {:error -0.022437763950071153, :complexity 24, :genotype [[0 1 0 0] [0 0 1 0] [2 3 0 5] [4 1 3 2] [0 0 0 1] [4 3 1 1]], :age 1} {:error 18.502165076507065, :complexity 8, :genotype [[0 0 1 0] [0 5 3 0]], :age 1} {:error -2.925055521072637E-4, :complexity 12, :genotype [[1 0 0 0] [0 0 1 0] [0 2 7 1]], :age 1} {:error 0.0, :complexity 16, :genotype [[1 0 0 0] [0 1 0 0] [0 0 1 0] [0 0 0 1]], :age 1} {:error 1.5832220019188459, :complexity 4, :genotype [[0 1 0 0]], :age 1} {:error 0.0, :complexity 8, :genotype [[1 0 0 0] [0 1 0 0]], :age 1} {:error 16.212236546842988, :complexity 4, :genotype [[0 0 1 0]], :age 1} {:error 1.5832220019188459, :complexity 4, :genotype [[1 0 0 0]], :age 1} {:error 0.0, :complexity 8, :genotype [[1 0 0 0] [0 1 0 0]], :age 1} {:error -0.03735029749918881, :complexity 16, :genotype [[1 0 0 0] [0 1 0 0] [0 0 1 0] [2 3 0 4]], :age 1} {:error -8.778546303851269E-4, :complexity 8, :genotype [[1 0 0 0] [2 3 0 3]], :age 1} {:error -4.813730779767095E-5, :complexity 12, :genotype [[1 0 0 0] [0 0 1 0] [0 1 2 5]], :age 1} {:error -0.0016882428922139103, :complexity 8, :genotype [[1 0 0 0] [5 1 1 1]], :age 1} {:error -0.030141422585557676, :complexity 12, :genotype [[1 0 0 0] [0 1 0 0] [0 3 4 0]], :age 1} {:error -1.0600768940669977E-5, :complexity 16, :genotype [[0 1 0 0] [0 0 1 0] [3 1 4 2] [0 0 0 1]], :age 1} {:error 0.028232322722737884, :complexity 12, :genotype [[1 0 0 0] [0 0 1 0] [0 0 0 1]], :age 1} {:error 0.0, :complexity 8, :genotype [[1 0 0 0] [0 1 0 0]], :age 1} {:error 1.5832220019188459, :complexity 4, :genotype [[1 0 0 0]], :age 1} {:error 0.00986814613767775, :complexity 8, :genotype [[0 1 0 0] [0 0 0 1]], :age 1} {:error -3.2655088756706024E-5, :complexity 16, :genotype [[2 4 2 1] [1 0 0 0] [0 1 0 0] [0 0 1 0]], :age 1} {:error -0.8745446366682642, :complexity 16, :genotype [[1 0 0 0] [6 0 0 0] [0 0 0 1] [4 0 1 1]], :age 1} {:error -0.006806603547829176, :complexity 16, :genotype [[1 0 0 0] [0 1 0 0] [7 0 0 3] [0 0 0 1]], :age 1} {:error 0.00986814613767775, :complexity 8, :genotype [[1 0 0 0] [0 0 1 0]], :age 1} {:error 1.5832220019188459, :complexity 4, :genotype [[0 1 0 0]], :age 1} {:error -5.058566372768392E-4, :complexity 24, :genotype [[1 0 0 0] [0 0 1 0] [5 0 3 0] [1 4 0 1] [0 0 0 1] [0 0 4 0]], :age 1} {:error -0.0014107972159827051, :complexity 12, :genotype [[0 1 0 0] [0 0 1 0] [6 0 2 2]], :age 1} {:error 0.00986814613767775, :complexity 8, :genotype [[1 0 0 0] [0 0 0 1]], :age 1} {:error 1.5832220019188459, :complexity 4, :genotype [[1 0 0 0]], :age 1} {:error 1.5832220019188459, :complexity 4, :genotype [[1 0 0 0]], :age 1} {:error -1.8146670833770209, :complexity 24, :genotype [[0 1 0 0] [0 0 1 0] [0 3 1 5] [2 5 2 1] [2 3 0 0] [0 0 0 1]], :age 1} {:error -5.024041690063155E-4, :complexity 16, :genotype [[1 0 0 0] [6 0 0 4] [2 0 0 7] [0 0 0 1]], :age 1} {:error 1.5832220019188459, :complexity 4, :genotype [[0 1 0 0]], :age 1} {:error -2.336597752106272E-4, :complexity 12, :genotype [[1 0 0 0] [0 0 0 1] [0 4 1 2]], :age 1} {:error 0.00986814613767775, :complexity 8, :genotype [[1 0 0 0] [0 0 1 0]], :age 1} {:error 0.0, :complexity 8, :genotype [[1 0 0 0] [0 1 0 0]], :age 1})
;;; ({:spea2-fitness 0.1587728278610361, :spea2-raw-fitness 0, :spea2-density 0.1587728278610361, :error 1.5832220019188459, :complexity 4, :genotype [[1 0 0 0]]} {:spea2-fitness 0.1587728278610361, :spea2-raw-fitness 0, :spea2-density 0.1587728278610361, :error 1.5832220019188459, :complexity 4, :genotype [[1 0 0 0]]} {:spea2-fitness 0.1587728278610361, :spea2-raw-fitness 0, :spea2-density 0.1587728278610361, :error 1.5832220019188459, :complexity 4, :genotype [[1 0 0 0]]} {:spea2-fitness 0.1587728278610361, :spea2-raw-fitness 0, :spea2-density 0.1587728278610361, :error 1.5832220019188459, :complexity 4, :genotype [[1 0 0 0]]} {:spea2-fitness 0.1587728278610361, :spea2-raw-fitness 0, :spea2-density 0.1587728278610361, :error 1.5832220019188459, :complexity 4, :genotype [[0 1 0 0]]} {:spea2-fitness 0.1587728278610361, :spea2-raw-fitness 0, :spea2-density 0.1587728278610361, :error 1.5832220019188459, :complexity 4, :genotype [[1 0 0 0]]} {:spea2-fitness 0.1843752345714079, :spea2-raw-fitness 0, :spea2-density 0.1843752345714079, :error -3.4282408092324927, :complexity 20, :genotype [[1 0 0 0] [0 1 0 0] [0 3 0 4] [0 0 1 0] [3 5 0 0]]} {:spea2-fitness 0.47474856856520636, :spea2-raw-fitness 0, :spea2-density 0.47474856856520636, :error -0.10715149932514655, :complexity 16, :genotype [[1 0 0 0] [0 0 1 0] [1 6 0 2] [0 1 0 0]]} {:spea2-fitness 0.497401033742391, :spea2-raw-fitness 0, :spea2-density 0.497401033742391, :error -0.010450184383633588, :complexity 12, :genotype [[0 3 6 0] [1 0 0 0] [0 1 0 0]]} {:spea2-fitness 0.49754507623879674, :spea2-raw-fitness 0, :spea2-density 0.49754507623879674, :error 0.0, :complexity 8, :genotype [[1 0 0 0] [0 1 0 0]]} {:spea2-fitness 0.49754507623879674, :spea2-raw-fitness 0, :spea2-density 0.49754507623879674, :error 0.0, :complexity 8, :genotype [[1 0 0 0] [0 1 0 0]]} {:spea2-fitness 0.49754507623879674, :spea2-raw-fitness 0, :spea2-density 0.49754507623879674, :error 0.0, :complexity 8, :genotype [[1 0 0 0] [0 1 0 0]]} {:spea2-fitness 0.49754507623879674, :spea2-raw-fitness 0, :spea2-density 0.49754507623879674, :error 0.0, :complexity 8, :genotype [[0 1 0 0] [1 0 0 0]]} {:spea2-fitness 0.49754507623879674, :spea2-raw-fitness 0, :spea2-density 0.49754507623879674, :error 0.0, :complexity 8, :genotype [[0 1 0 0] [1 0 0 0]]} {:spea2-fitness 0.49754507623879674, :spea2-raw-fitness 0, :spea2-density 0.49754507623879674, :error 0.0, :complexity 8, :genotype [[1 0 0 0] [0 1 0 0]]} {:spea2-fitness 23.156573566287626, :spea2-raw-fitness 23, :spea2-density 0.15657356628762706, :error -1.8146670833770209, :complexity 24, :genotype [[0 0 0 1] [0 0 1 0] [2 5 2 1] [0 1 0 0] [2 3 0 0] [0 3 1 5]]} {:spea2-fitness 23.266181197165032, :spea2-raw-fitness 23, :spea2-density 0.26618119716503075, :error -1.7613582620178936, :complexity 20, :genotype [[1 3 3 3] [6 1 0 0] [0 0 1 0] [0 1 0 0] [0 0 0 1]]} {:spea2-fitness 43.48217309378113, :spea2-raw-fitness 43, :spea2-density 0.4821730937811286, :error -0.0747173844611858, :complexity 16, :genotype [[0 0 1 0] [0 0 0 1] [1 6 0 2] [0 1 0 0]]} {:spea2-fitness 44.29101699904094, :spea2-raw-fitness 44, :spea2-density 0.29101699904094036, :error -1.4399696055586102, :complexity 20, :genotype [[0 0 0 1] [0 0 1 0] [0 1 0 0] [1 0 0 0] [3 1 0 0]]} {:spea2-fitness 58.49885669839754, :spea2-raw-fitness 58, :spea2-density 0.49885669839753904, :error -0.004583687484335895, :complexity 12, :genotype [[2 3 0 3] [1 0 0 0] [0 1 0 0]]} {:spea2-fitness 69.16410998142895, :spea2-raw-fitness 69, :spea2-density 0.16410998142895125, :error -0.8765963074677557, :complexity 24, :genotype [[0 0 0 1] [0 0 1 0] [4 0 1 1] [1 0 0 0] [2 4 2 1] [6 0 0 0]]} {:spea2-fitness 85.49084472790986, :spea2-raw-fitness 85, :spea2-density 0.4908447279098561, :error -0.03735029749918881, :complexity 16, :genotype [[0 1 0 0] [1 0 0 0] [2 3 0 4] [0 0 1 0]]} {:spea2-fitness 90.1513066067542, :spea2-raw-fitness 90, :spea2-density 0.15130660675420074, :error 16.212236546842988, :complexity 4, :genotype [[0 0 1 0]]} {:spea2-fitness 90.1513066067542, :spea2-raw-fitness 90, :spea2-density 0.15130660675420074, :error 16.212236546842988, :complexity 4, :genotype [[0 0 1 0]]} {:spea2-fitness 90.1513066067542, :spea2-raw-fitness 90, :spea2-density 0.15130660675420074, :error 16.212236546842988, :complexity 4, :genotype [[0 0 1 0]]} {:spea2-fitness 90.1513066067542, :spea2-raw-fitness 90, :spea2-density 0.15130660675420074, :error 16.212236546842988, :complexity 4, :genotype [[0 0 1 0]]} {:spea2-fitness 90.1513066067542, :spea2-raw-fitness 90, :spea2-density 0.15130660675420074, :error 16.212236546842988, :complexity 4, :genotype [[0 0 0 1]]} {:spea2-fitness 90.1513066067542, :spea2-raw-fitness 90, :spea2-density 0.15130660675420074, :error 16.212236546842988, :complexity 4, :genotype [[0 0 1 0]]} {:spea2-fitness 90.1513066067542, :spea2-raw-fitness 90, :spea2-density 0.15130660675420074, :error 16.212236546842988, :complexity 4, :genotype [[0 0 1 0]]} {:spea2-fitness 112.49959208486791, :spea2-raw-fitness 112, :spea2-density 0.49959208486791873, :error -0.0016329927732507256, :complexity 12, :genotype [[1 0 0 0] [2 0 0 1] [4 2 0 2]]} {:spea2-fitness 125.49258788316196, :spea2-raw-fitness 125, :spea2-density 0.49258788316196694, :error -0.03014074658160298, :complexity 16, :genotype [[0 1 0 0] [0 0 0 1] [0 3 4 0] [0 0 1 0]]} {:spea2-fitness 149.48847066590756, :spea2-raw-fitness 149, :spea2-density 0.48847066590756033, :error -0.04899291061013166, :complexity 20, :genotype [[0 0 1 0] [1 0 0 0] [0 0 0 1] [1 0 6 0] [2 5 2 0]]} {:spea2-fitness 157.4998411652939, :spea2-raw-fitness 157, :spea2-density 0.4998411652939107, :error -6.355407162029653E-4, :complexity 12, :genotype [[1 0 0 0] [1 0 7 2] [0 1 0 0]]} {:spea2-fitness 164.4968906450108, :spea2-raw-fitness 164, :spea2-density 0.49689064501080016, :error -0.012524971127530439, :complexity 16, :genotype [[0 1 0 0] [1 0 0 0] [1 0 0 6] [0 0 0 1]]} {:spea2-fitness 197.4999142559156, :spea2-raw-fitness 197, :spea2-density 0.49991425591559413, :error -2.925055521072637E-4, :complexity 12, :genotype [[0 2 7 1] [0 0 1 0] [1 0 0 0]]} {:spea2-fitness 222.4975450762388, :spea2-raw-fitness 222, :spea2-density 0.49754507623879674, :error 0.00986814613767775, :complexity 8, :genotype [[0 0 0 1] [0 1 0 0]]} {:spea2-fitness 222.4975450762388, :spea2-raw-fitness 222, :spea2-density 0.49754507623879674, :error 0.00986814613767775, :complexity 8, :genotype [[1 0 0 0] [0 0 1 0]]} {:spea2-fitness 222.4975450762388, :spea2-raw-fitness 222, :spea2-density 0.49754507623879674, :error 0.00986814613767775, :complexity 8, :genotype [[1 0 0 0] [0 0 1 0]]} {:spea2-fitness 222.4975450762388, :spea2-raw-fitness 222, :spea2-density 0.49754507623879674, :error 0.00986814613767775, :complexity 8, :genotype [[1 0 0 0] [0 0 0 1]]} {:spea2-fitness 231.4998732694335, :spea2-raw-fitness 231, :spea2-density 0.49987326943348304, :error -1.2848993246911792E-4, :complexity 12, :genotype [[0 1 0 0] [1 0 0 0] [3 1 4 2]]} {:spea2-fitness 246.4969755158226, :spea2-raw-fitness 246, :spea2-density 0.49697551582259897, :error -0.013615412403065573, :complexity 20, :genotype [[0 0 0 1] [0 1 0 0] [1 0 0 0] [0 0 4 3] [0 0 1 0]]} {:spea2-fitness 255.09999980807544, :spea2-raw-fitness 255, :spea2-density 0.09999980807542892, :error -0.023512378913830134, :complexity 28, :genotype [[2 3 3 2] [2 3 0 5] [0 0 0 1] [0 1 7 0] [4 3 1 1] [0 1 0 0] [2 6 0 2]]} {:spea2-fitness 257.0999999317823, :spea2-raw-fitness 257, :spea2-density 0.09999993178230955, :error -0.016436085323039997, :complexity 28, :genotype [[1 4 0 1] [0 0 1 0] [5 0 3 0] [0 0 4 0] [1 0 0 0] [0 0 0 1] [0 1 0 0]]} {:spea2-fitness 270.166666410893, :spea2-raw-fitness 270, :spea2-density 0.16666641089302964, :error -0.01310161150949659, :complexity 24, :genotype [[0 0 0 1] [1 0 0 0] [0 1 0 0] [0 5 4 1] [0 0 1 0] [1 2 0 6]]} {:spea2-fitness 310.4989476605517, :spea2-raw-fitness 310, :spea2-density 0.498947660551744, :error -0.004218235824945162, :complexity 16, :genotype [[1 0 3 0] [0 1 0 0] [1 0 0 0] [1 3 0 2]]} {:spea2-fitness 338.4989696452205, :spea2-raw-fitness 338, :spea2-density 0.49896964522050846, :error -0.004129929703584058, :complexity 16, :genotype [[0 0 1 0] [1 0 0 0] [1 1 3 4] [5 2 1 0]]} {:spea2-fitness 353.49780063711376, :spea2-raw-fitness 353, :spea2-density 0.4978006371137671, :error -0.010280170393330347, :complexity 20, :genotype [[0 0 0 1] [0 0 1 0] [1 0 5 4] [1 0 0 0] [4 0 0 6]]} {:spea2-fitness 365.4995646938797, :spea2-raw-fitness 365, :spea2-density 0.49956469387968977, :error -0.0017427417335262538, :complexity 16, :genotype [[5 1 1 1] [0 0 1 0] [0 1 2 5] [1 0 0 0]]} {:spea2-fitness 366.49858046730265, :spea2-raw-fitness 366, :spea2-density 0.4985804673026569, :error -0.006806603547829176, :complexity 20, :genotype [[0 0 1 0] [7 0 0 3] [1 0 0 0] [0 0 0 1] [0 1 0 0]]} {:spea2-fitness 378.4987838734761, :spea2-raw-fitness 378, :spea2-density 0.49878387347612557, :error -0.005988672929277489, :complexity 20, :genotype [[0 0 0 1] [1 0 0 0] [6 0 0 4] [2 0 0 7] [0 1 0 0]]} {:spea2-fitness 418.1666666666667, :spea2-raw-fitness 418, :spea2-density 0.16666666666666666, :error 16.212236546842988, :complexity 8, :genotype [[0 0 1 0] [0 0 0 1]]} {:spea2-fitness 423.16578719678836, :spea2-raw-fitness 423, :spea2-density 0.165787196788352, :error 16.71784908976807, :complexity 8, :genotype [[0 0 0 1] [7 0 2 1]]} {:spea2-fitness 427.1513066067542, :spea2-raw-fitness 427, :spea2-density 0.15130660675420074, :error 18.502165076507065, :complexity 8, :genotype [[0 0 1 0] [0 5 3 0]]} {:spea2-fitness 433.4996579060447, :spea2-raw-fitness 433, :spea2-density 0.49965790604474536, :error -0.001369312688205371, :complexity 16, :genotype [[0 0 0 1] [1 0 0 0] [0 2 6 1] [0 0 1 0]]} {:spea2-fitness 443.4990149977923, :spea2-raw-fitness 443, :spea2-density 0.4990149977922802, :error -0.004564171439618814, :complexity 20, :genotype [[6 0 2 2] [0 1 0 0] [1 0 0 0] [0 0 1 0] [0 0 0 1]]} {:spea2-fitness 453.4990262723513, :spea2-raw-fitness 453, :spea2-density 0.4990262723512904, :error -0.004518896012615737, :complexity 20, :genotype [[1 0 0 0] [0 0 1 0] [0 1 0 0] [0 0 0 1] [0 5 3 0]]} {:spea2-fitness 454.49980673131074, :spea2-raw-fitness 454, :spea2-density 0.4998067313107236, :error -7.733736949462403E-4, :complexity 16, :genotype [[3 0 0 7] [1 0 0 0] [0 0 1 0] [0 1 0 0]]} {:spea2-fitness 486.4998411652939, :spea2-raw-fitness 486, :spea2-density 0.4998411652939107, :error 0.0, :complexity 12, :genotype [[0 0 1 0] [1 0 0 0] [0 1 0 0]]} {:spea2-fitness 486.4998411652939, :spea2-raw-fitness 486, :spea2-density 0.4998411652939107, :error 0.0, :complexity 12, :genotype [[1 0 0 0] [0 0 1 0] [0 1 0 0]]} {:spea2-fitness 486.4998411652939, :spea2-raw-fitness 486, :spea2-density 0.4998411652939107, :error 0.0, :complexity 12, :genotype [[0 1 0 0] [3 0 4 0] [0 0 1 0]]} {:spea2-fitness 486.4998411652939, :spea2-raw-fitness 486, :spea2-density 0.4998411652939107, :error 0.0, :complexity 12, :genotype [[0 1 0 0] [0 0 1 0] [1 0 0 0]]} {:spea2-fitness 486.4998411652939, :spea2-raw-fitness 486, :spea2-density 0.4998411652939107, :error 0.0, :complexity 12, :genotype [[1 0 0 0] [0 0 0 1] [0 1 0 0]]} {:spea2-fitness 486.4998411652939, :spea2-raw-fitness 486, :spea2-density 0.4998411652939107, :error 0.0, :complexity 12, :genotype [[0 1 0 0] [1 0 0 0] [0 0 1 0]]} {:spea2-fitness 486.4998411652939, :spea2-raw-fitness 486, :spea2-density 0.4998411652939107, :error 0.0, :complexity 12, :genotype [[0 0 1 0] [0 1 0 0] [1 0 0 0]]} {:spea2-fitness 513.4998696019479, :spea2-raw-fitness 513, :spea2-density 0.499869601947911, :error -5.217282730569098E-4, :complexity 16, :genotype [[0 4 1 2] [1 0 0 0] [0 0 0 1] [0 1 0 0]]} {:spea2-fitness 517.4992192626005, :spea2-raw-fitness 517, :spea2-density 0.49921926260045923, :error -0.003744219084468871, :complexity 20, :genotype [[2 4 0 1] [0 0 1 0] [1 0 0 0] [0 1 0 0] [0 0 0 1]]} {:spea2-fitness 525.4994337954942, :spea2-raw-fitness 525, :spea2-density 0.49943379549410116, :error -0.0025642878594982982, :complexity 20, :genotype [[0 0 0 1] [0 2 6 0] [0 0 1 0] [0 1 0 0] [1 1 4 3]]} {:spea2-fitness 532.4993066874484, :spea2-raw-fitness 532, :spea2-density 0.49930668744843565, :error -0.0017870704353952817, :complexity 20, :genotype [[0 1 0 0] [0 0 0 1] [0 0 1 0] [1 0 0 0] [4 3 3 0]]} {:spea2-fitness 597.4999884627923, :spea2-raw-fitness 597, :spea2-density 0.49998846279226716, :error -4.614989581352499E-5, :complexity 16, :genotype [[1 0 0 0] [0 1 0 0] [1 3 3 3] [0 0 0 1]]} {:spea2-fitness 606.4992211348696, :spea2-raw-fitness 606, :spea2-density 0.4992211348696282, :error -0.0014438502994683718, :complexity 20, :genotype [[0 0 1 0] [0 0 0 1] [4 0 2 4] [1 0 0 0] [6 0 4 0]]} {:spea2-fitness 610.4999897260692, :spea2-raw-fitness 610, :spea2-density 0.49998972606914766, :error -4.1096567855760434E-5, :complexity 16, :genotype [[0 2 1 7] [0 0 1 0] [0 1 0 0] [0 0 0 1]]} {:spea2-fitness 622.4999975693945, :spea2-raw-fitness 622, :spea2-density 0.49999756939456175, :error -9.722469016059903E-6, :complexity 16, :genotype [[0 0 0 1] [0 1 0 0] [0 0 1 0] [1 5 2 0]]} {:spea2-fitness 632.4991385205645, :spea2-raw-fitness 632, :spea2-density 0.49913852056455293, :error -0.0011123062761258483, :complexity 20, :genotype [[0 0 1 0] [0 1 0 0] [1 1 4 2] [1 0 0 0] [0 0 0 1]]} {:spea2-fitness 695.4990149977923, :spea2-raw-fitness 695, :spea2-density 0.4990149977922802, :error -6.163854529145253E-4, :complexity 20, :genotype [[0 0 0 1] [0 5 1 1] [1 0 0 0] [2 0 6 0] [2 0 4 2]]} {:spea2-fitness 711.4930401654666, :spea2-raw-fitness 711, :spea2-density 0.4930401654666369, :error 0.028232322722737884, :complexity 12, :genotype [[0 0 0 1] [0 0 1 0] [0 1 0 0]]} {:spea2-fitness 711.4930401654666, :spea2-raw-fitness 711, :spea2-density 0.4930401654666369, :error 0.028232322722737884, :complexity 12, :genotype [[0 0 1 0] [0 1 0 0] [0 0 0 1]]} {:spea2-fitness 711.4930401654666, :spea2-raw-fitness 711, :spea2-density 0.4930401654666369, :error 0.028232322722737884, :complexity 12, :genotype [[0 0 1 0] [0 1 0 0] [0 0 0 1]]} {:spea2-fitness 711.4930401654666, :spea2-raw-fitness 711, :spea2-density 0.4930401654666369, :error 0.028232322722737884, :complexity 12, :genotype [[0 0 1 0] [0 0 0 1] [0 1 0 0]]} {:spea2-fitness 711.4930401654666, :spea2-raw-fitness 711, :spea2-density 0.4930401654666369, :error 0.028232322722737884, :complexity 12, :genotype [[0 0 1 0] [0 0 0 1] [1 0 0 0]]} {:spea2-fitness 711.4930401654666, :spea2-raw-fitness 711, :spea2-density 0.4930401654666369, :error 0.028232322722737884, :complexity 12, :genotype [[1 0 0 0] [0 0 1 0] [0 0 0 1]]} {:spea2-fitness 711.4930401654666, :spea2-raw-fitness 711, :spea2-density 0.4930401654666369, :error 0.028232322722737884, :complexity 12, :genotype [[0 0 1 0] [0 1 0 0] [0 0 0 1]]} {:spea2-fitness 711.4930401654666, :spea2-raw-fitness 711, :spea2-density 0.4930401654666369, :error 0.028232322722737884, :complexity 12, :genotype [[0 0 1 0] [0 0 0 1] [0 1 0 0]]} {:spea2-fitness 715.4989699756727, :spea2-raw-fitness 715, :spea2-density 0.4989699756727261, :error -4.355690086403119E-4, :complexity 20, :genotype [[1 0 0 0] [0 4 2 4] [1 5 1 3] [0 0 1 0] [0 2 2 4]]} {:spea2-fitness 717.4989354540492, :spea2-raw-fitness 717, :spea2-density 0.4989354540491706, :error -2.9690222798052447E-4, :complexity 20, :genotype [[0 1 0 0] [1 0 3 1] [0 6 1 0] [3 3 3 1] [0 0 0 1]]} {:spea2-fitness 730.0999999996911, :spea2-raw-fitness 730, :spea2-density 0.09999999969117031, :error -4.093642274742216E-4, :complexity 28, :genotype [[0 1 0 0] [3 3 3 1] [0 0 2 2] [0 0 1 0] [1 1 6 2] [1 0 0 0] [0 0 0 1]]} {:spea2-fitness 841.1665174384513, :spea2-raw-fitness 841, :spea2-density 0.1665174384512523, :error 4.005377030182627, :complexity 12, :genotype [[1 5 2 0] [3 5 0 1] [0 0 0 1]]} {:spea2-fitness 845.0888760523067, :spea2-raw-fitness 845, :spea2-density 0.08887605230660797, :error 11.565456733340245, :complexity 12, :genotype [[0 0 1 0] [0 4 2 4] [1 5 0 2]]} {:spea2-fitness 902.0608963008054, :spea2-raw-fitness 902, :spea2-density 0.060896300805457056, :error 28.211219304298957, :complexity 12, :genotype [[0 0 1 0] [0 0 0 1] [4 1 3 2]]} {:spea2-fitness 904.0468778789606, :spea2-raw-fitness 904, :spea2-density 0.04687787896060714, :error 33.81130242328212, :complexity 12, :genotype [[0 0 0 1] [0 0 1 0] [2 0 5 0]]} {:spea2-fitness 1016.5, :spea2-raw-fitness 1016, :spea2-density 0.5, :error 0.0, :complexity 16, :genotype [[0 0 1 0] [0 1 0 0] [1 0 0 0] [0 0 0 1]]} {:spea2-fitness 1016.5, :spea2-raw-fitness 1016, :spea2-density 0.5, :error 0.0, :complexity 16, :genotype [[0 0 0 1] [0 0 1 0] [1 0 0 0] [0 1 0 0]]} {:spea2-fitness 1016.5, :spea2-raw-fitness 1016, :spea2-density 0.5, :error 0.0, :complexity 16, :genotype [[1 0 0 0] [0 1 0 0] [0 0 0 1] [0 0 1 0]]} {:spea2-fitness 1016.5, :spea2-raw-fitness 1016, :spea2-density 0.5, :error 0.0, :complexity 16, :genotype [[1 0 0 0] [0 0 1 0] [0 0 0 1] [0 1 0 0]]} {:spea2-fitness 1016.5, :spea2-raw-fitness 1016, :spea2-density 0.5, :error 0.0, :complexity 16, :genotype [[0 1 0 0] [0 0 0 1] [0 0 1 0] [1 0 0 0]]} {:spea2-fitness 1016.5, :spea2-raw-fitness 1016, :spea2-density 0.5, :error 0.0, :complexity 16, :genotype [[1 0 0 0] [0 0 0 1] [0 0 1 0] [0 1 0 0]]} {:spea2-fitness 1016.5, :spea2-raw-fitness 1016, :spea2-density 0.5, :error 0.0, :complexity 16, :genotype [[1 0 0 0] [0 0 1 0] [0 0 0 1] [0 1 0 0]]} {:spea2-fitness 1016.5, :spea2-raw-fitness 1016, :spea2-density 0.5, :error 0.0, :complexity 16, :genotype [[0 0 0 1] [1 0 0 0] [0 1 0 0] [0 0 1 0]]} {:spea2-fitness 1016.5, :spea2-raw-fitness 1016, :spea2-density 0.5, :error 0.0, :complexity 16, :genotype [[0 0 1 0] [0 0 0 1] [0 1 0 0] [1 0 0 0]]} {:spea2-fitness 1016.5, :spea2-raw-fitness 1016, :spea2-density 0.5, :error 0.0, :complexity 16, :genotype [[1 0 0 0] [0 1 0 0] [0 0 1 0] [0 0 0 1]]} {:spea2-fitness 1284.0374896040128, :spea2-raw-fitness 1284, :spea2-density 0.037489604012653506, :error 37.771673155871085, :complexity 16, :genotype [[1 1 4 2] [1 0 8 0] [0 0 0 1] [0 0 1 0]]})
;;; .&quot;Elapsed time: 8128.367312 msecs&quot;
;;; 
;; <-
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/result</span>","value":"#'user/result"}
;; <=

;; @@
(.DestroyMathKernel client)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>true</span>","value":"true"}
;; <=

;; @@
(.DisconnectServer user/client)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@
(mapv #(println (:genotype %)) (sort-by :error (:elite result)))
;; @@
;; ->
;;; [[1 0 0 0] [0 1 0 0] [0 3 0 4] [0 0 1 0] [3 5 0 0]]
;;; [[0 0 0 1] [0 0 1 0] [2 5 2 1] [0 1 0 0] [2 3 0 0] [0 3 1 5]]
;;; [[1 3 3 3] [6 1 0 0] [0 0 1 0] [0 1 0 0] [0 0 0 1]]
;;; [[0 0 0 1] [0 0 1 0] [0 1 0 0] [1 0 0 0] [3 1 0 0]]
;;; [[0 0 0 1] [0 0 1 0] [4 0 1 1] [1 0 0 0] [2 4 2 1] [6 0 0 0]]
;;; [[1 0 0 0] [0 0 1 0] [1 6 0 2] [0 1 0 0]]
;;; [[0 0 1 0] [0 0 0 1] [1 6 0 2] [0 1 0 0]]
;;; [[0 0 1 0] [1 0 0 0] [0 0 0 1] [1 0 6 0] [2 5 2 0]]
;;; [[0 1 0 0] [1 0 0 0] [2 3 0 4] [0 0 1 0]]
;;; [[0 1 0 0] [0 0 0 1] [0 3 4 0] [0 0 1 0]]
;;; [[2 3 3 2] [2 3 0 5] [0 0 0 1] [0 1 7 0] [4 3 1 1] [0 1 0 0] [2 6 0 2]]
;;; [[1 4 0 1] [0 0 1 0] [5 0 3 0] [0 0 4 0] [1 0 0 0] [0 0 0 1] [0 1 0 0]]
;;; [[0 0 0 1] [0 1 0 0] [1 0 0 0] [0 0 4 3] [0 0 1 0]]
;;; [[0 0 0 1] [1 0 0 0] [0 1 0 0] [0 5 4 1] [0 0 1 0] [1 2 0 6]]
;;; [[0 1 0 0] [1 0 0 0] [1 0 0 6] [0 0 0 1]]
;;; [[0 3 6 0] [1 0 0 0] [0 1 0 0]]
;;; [[0 0 0 1] [0 0 1 0] [1 0 5 4] [1 0 0 0] [4 0 0 6]]
;;; [[0 0 1 0] [7 0 0 3] [1 0 0 0] [0 0 0 1] [0 1 0 0]]
;;; [[0 0 0 1] [1 0 0 0] [6 0 0 4] [2 0 0 7] [0 1 0 0]]
;;; [[2 3 0 3] [1 0 0 0] [0 1 0 0]]
;;; [[6 0 2 2] [0 1 0 0] [1 0 0 0] [0 0 1 0] [0 0 0 1]]
;;; [[1 0 0 0] [0 0 1 0] [0 1 0 0] [0 0 0 1] [0 5 3 0]]
;;; [[1 0 3 0] [0 1 0 0] [1 0 0 0] [1 3 0 2]]
;;; [[0 0 1 0] [1 0 0 0] [1 1 3 4] [5 2 1 0]]
;;; [[2 4 0 1] [0 0 1 0] [1 0 0 0] [0 1 0 0] [0 0 0 1]]
;;; [[0 0 0 1] [0 2 6 0] [0 0 1 0] [0 1 0 0] [1 1 4 3]]
;;; [[0 1 0 0] [0 0 0 1] [0 0 1 0] [1 0 0 0] [4 3 3 0]]
;;; [[5 1 1 1] [0 0 1 0] [0 1 2 5] [1 0 0 0]]
;;; [[1 0 0 0] [2 0 0 1] [4 2 0 2]]
;;; [[0 0 1 0] [0 0 0 1] [4 0 2 4] [1 0 0 0] [6 0 4 0]]
;;; [[0 0 0 1] [1 0 0 0] [0 2 6 1] [0 0 1 0]]
;;; [[0 0 1 0] [0 1 0 0] [1 1 4 2] [1 0 0 0] [0 0 0 1]]
;;; [[3 0 0 7] [1 0 0 0] [0 0 1 0] [0 1 0 0]]
;;; [[1 0 0 0] [1 0 7 2] [0 1 0 0]]
;;; [[0 0 0 1] [0 5 1 1] [1 0 0 0] [2 0 6 0] [2 0 4 2]]
;;; [[0 4 1 2] [1 0 0 0] [0 0 0 1] [0 1 0 0]]
;;; [[1 0 0 0] [0 4 2 4] [1 5 1 3] [0 0 1 0] [0 2 2 4]]
;;; [[0 1 0 0] [3 3 3 1] [0 0 2 2] [0 0 1 0] [1 1 6 2] [1 0 0 0] [0 0 0 1]]
;;; [[0 1 0 0] [1 0 3 1] [0 6 1 0] [3 3 3 1] [0 0 0 1]]
;;; [[0 2 7 1] [0 0 1 0] [1 0 0 0]]
;;; [[0 1 0 0] [1 0 0 0] [3 1 4 2]]
;;; [[1 0 0 0] [0 1 0 0] [1 3 3 3] [0 0 0 1]]
;;; [[0 2 1 7] [0 0 1 0] [0 1 0 0] [0 0 0 1]]
;;; [[0 0 0 1] [0 1 0 0] [0 0 1 0] [1 5 2 0]]
;;; [[1 0 0 0] [0 1 0 0]]
;;; [[1 0 0 0] [0 1 0 0]]
;;; [[1 0 0 0] [0 1 0 0]]
;;; [[0 1 0 0] [1 0 0 0]]
;;; [[0 1 0 0] [1 0 0 0]]
;;; [[1 0 0 0] [0 1 0 0]]
;;; [[0 0 1 0] [1 0 0 0] [0 1 0 0]]
;;; [[1 0 0 0] [0 0 1 0] [0 1 0 0]]
;;; [[0 1 0 0] [3 0 4 0] [0 0 1 0]]
;;; [[0 1 0 0] [0 0 1 0] [1 0 0 0]]
;;; [[1 0 0 0] [0 0 0 1] [0 1 0 0]]
;;; [[0 1 0 0] [1 0 0 0] [0 0 1 0]]
;;; [[0 0 1 0] [0 1 0 0] [1 0 0 0]]
;;; [[0 0 1 0] [0 1 0 0] [1 0 0 0] [0 0 0 1]]
;;; [[0 0 0 1] [0 0 1 0] [1 0 0 0] [0 1 0 0]]
;;; [[1 0 0 0] [0 1 0 0] [0 0 0 1] [0 0 1 0]]
;;; [[1 0 0 0] [0 0 1 0] [0 0 0 1] [0 1 0 0]]
;;; [[0 1 0 0] [0 0 0 1] [0 0 1 0] [1 0 0 0]]
;;; [[1 0 0 0] [0 0 0 1] [0 0 1 0] [0 1 0 0]]
;;; [[1 0 0 0] [0 0 1 0] [0 0 0 1] [0 1 0 0]]
;;; [[0 0 0 1] [1 0 0 0] [0 1 0 0] [0 0 1 0]]
;;; [[0 0 1 0] [0 0 0 1] [0 1 0 0] [1 0 0 0]]
;;; [[1 0 0 0] [0 1 0 0] [0 0 1 0] [0 0 0 1]]
;;; [[0 0 0 1] [0 1 0 0]]
;;; [[1 0 0 0] [0 0 1 0]]
;;; [[1 0 0 0] [0 0 1 0]]
;;; [[1 0 0 0] [0 0 0 1]]
;;; [[0 0 0 1] [0 0 1 0] [0 1 0 0]]
;;; [[0 0 1 0] [0 1 0 0] [0 0 0 1]]
;;; [[0 0 1 0] [0 1 0 0] [0 0 0 1]]
;;; [[0 0 1 0] [0 0 0 1] [0 1 0 0]]
;;; [[0 0 1 0] [0 0 0 1] [1 0 0 0]]
;;; [[1 0 0 0] [0 0 1 0] [0 0 0 1]]
;;; [[0 0 1 0] [0 1 0 0] [0 0 0 1]]
;;; [[0 0 1 0] [0 0 0 1] [0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; [[1 0 0 0]]
;;; [[0 1 0 0]]
;;; [[1 0 0 0]]
;;; [[1 5 2 0] [3 5 0 1] [0 0 0 1]]
;;; [[0 0 1 0] [0 4 2 4] [1 5 0 2]]
;;; [[0 0 1 0]]
;;; [[0 0 1 0]]
;;; [[0 0 1 0]]
;;; [[0 0 1 0]]
;;; [[0 0 0 1]]
;;; [[0 0 1 0]]
;;; [[0 0 1 0]]
;;; [[0 0 1 0] [0 0 0 1]]
;;; [[0 0 0 1] [7 0 2 1]]
;;; [[0 0 1 0] [0 5 3 0]]
;;; [[0 0 1 0] [0 0 0 1] [4 1 3 2]]
;;; [[0 0 0 1] [0 0 1 0] [2 0 5 0]]
;;; [[1 1 4 2] [1 0 8 0] [0 0 0 1] [0 0 1 0]]
;;; 
;; <-
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}],"value":"[nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil]"}
;; <=

;; @@

;; @@
