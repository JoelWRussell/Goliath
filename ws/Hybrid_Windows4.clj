;; gorilla-repl.fileformat = 1

;; **
;;; #HYBRID WITH VARIABLE DEGREES OF FREEDOM
;;; 
;;; The mma code has all been moved into LagrangeSolverCompact.m, the Java class is LagrangianScore.java. I have followed Jony's idea of memoizing the score function. This dramatically improves speed. What this does is that it keeps a record of previous calls to the score function and the corresponding result. So if it calls with the same arguments etc then it is able to find the answer via some sort of look up table.
;;; 
;;; It also saves each zeitgeist to disk - every generation.
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
            [criterium.core :as criterium]
            [semantic-csv.core :as csv])
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
(def OS (System/getProperty "os.name"))
(if (= OS "Mac OS X")
(def mathKernelSz "/Applications/Mathematica.app/Contents/MacOS/MathKernel")  
(def mathKernelSz "c:/program files/wolfram research/mathematica/10.0/mathkernel.exe"))

OS 
mathKernelSz
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-string'>&quot;/Applications/Mathematica.app/Contents/MacOS/MathKernel&quot;</span>","value":"\"/Applications/Mathematica.app/Contents/MacOS/MathKernel\""}
;; <=

;; **
;;; CHOOSE THE DATA - NEEDS TO BE IN RESOURCES/
;; **

;; @@
(def experimentalDataSz "sho_coupled_0.1_sim.csv")
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/experimentalDataSz</span>","value":"#'goliath/experimentalDataSz"}
;; <=

;; @@
(def dat (mapv #(vector (first %))
               (csv/slurp-csv "resources/sho_coupled_0.1_sim.csv"
  							  :header false
                              :cast-fns {0 #(Double/parseDouble %) 1 #(Double/parseDouble %)})))

(def dat2 (mapv #(vector (second %))
               (csv/slurp-csv "resources/sho_coupled_0.1_sim.csv"
                              :header false
                              :cast-fns {0 #(Double/parseDouble %) 1 #(Double/parseDouble %)})))

"Theta"

(plot/list-plot (mapv first dat) :joined true)

"Omega"

(plot/list-plot (mapv first dat2) :joined true)


;; @@
;; =>
;;; {"type":"vega","content":{"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"5ad73c2c-07da-49f8-bd98-1e4fc6917fa2","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":{"data":"5ad73c2c-07da-49f8-bd98-1e4fc6917fa2","field":"data.y"}}],"marks":[{"type":"line","from":{"data":"5ad73c2c-07da-49f8-bd98-1e4fc6917fa2"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"#FF29D2"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}}],"data":[{"name":"5ad73c2c-07da-49f8-bd98-1e4fc6917fa2","values":[{"x":0,"y":2.0},{"x":1,"y":1.9284804760741014},{"x":2,"y":1.847869707654912},{"x":3,"y":1.759147940170956},{"x":4,"y":1.6633539300874196},{"x":5,"y":1.5615627444976496},{"x":6,"y":1.4548630903699578},{"x":7,"y":1.3443348754972466},{"x":8,"y":1.2310275172785075},{"x":9,"y":1.1159395524551143},{"x":10,"y":1.0000000198584484},{"x":11,"y":0.884052161868249},{"x":12,"y":0.7688397122512814},{"x":13,"y":0.6549961985048567},{"x":14,"y":0.5430374695353829},{"x":15,"y":0.4333576046967714},{"x":16,"y":0.32622828094447165},{"x":17,"y":0.2218015812237824},{"x":18,"y":0.12011612593437726},{"x":19,"y":0.02110633326264301},{"x":20,"y":-0.07538547415089168},{"x":21,"y":-0.16959447211977985},{"x":22,"y":-0.2618166633840412},{"x":23,"y":-0.35239062947183336},{"x":24,"y":-0.4416781845088593},{"x":25,"y":-0.5300445077444007},{"x":26,"y":-0.6178383035773294},{"x":27,"y":-0.705372547526059},{"x":28,"y":-0.7929063568341025},{"x":29,"y":-0.88062849998643},{"x":30,"y":-0.9686430148760465},{"x":31,"y":-1.0569573420664895},{"x":32,"y":-1.1454733229161638},{"x":33,"y":-1.2339813292889525},{"x":34,"y":-1.3221577103472784},{"x":35,"y":-1.4095656548744435},{"x":36,"y":-1.4956594741653633},{"x":37,"y":-1.5797922236920523},{"x":38,"y":-1.6612264815843205},{"x":39,"y":-1.7391480398024644},{"x":40,"y":-1.8126821585831931},{"x":41,"y":-1.8809119854688685},{"x":42,"y":-1.9428986720936203},{"x":43,"y":-1.9977026777887568},{"x":44,"y":-2.044405720677371},{"x":45,"y":-2.082132804456243},{"x":46,"y":-2.110073758399698},{"x":47,"y":-2.1275037306372786},{"x":48,"y":-2.133802107707231},{"x":49,"y":-2.128469353974073},{"x":50,"y":-2.1111413434194244},{"x":51,"y":-2.0816007929888336},{"x":52,"y":-2.0397854974273417},{"x":53,"y":-1.9857931381250578},{"x":54,"y":-1.9198825296991564},{"x":55,"y":-1.8424712574666764},{"x":56,"y":-1.7541297447766149},{"x":57,"y":-1.655571898771775},{"x":58,"y":-1.5476425469860653},{"x":59,"y":-1.431301972250347},{"x":60,"y":-1.3076079261555191},{"x":61,"y":-1.1776955628380277},{"x":62,"y":-1.0427557895776662},{"x":63,"y":-0.9040125671716012},{"x":64,"y":-0.7626997233384684},{"x":65,"y":-0.6200378489656382},{"x":66,"y":-0.4772118474306306},{"x":67,"y":-0.33534968782029867},{"x":68,"y":-0.19550287980466313},{"x":69,"y":-0.05862913648959707},{"x":70,"y":0.07442234839953822},{"x":71,"y":0.2029226740930839},{"x":72,"y":0.3262717872830035},{"x":73,"y":0.44400362591251114},{"x":74,"y":0.5557877774505213},{"x":75,"y":0.6614276207541374},{"x":76,"y":0.760855066199606},{"x":77,"y":0.8541220585200456},{"x":78,"y":0.9413890972582583},{"x":79,"y":1.0229111375410067},{"x":80,"y":1.0990212412755835},{"x":81,"y":1.170112478892747},{"x":82,"y":1.2366185681636452},{"x":83,"y":1.2989937911787741},{"x":84,"y":1.3576927683000244},{"x":85,"y":1.4131506108921674},{"x":86,"y":1.4657640511313044},{"x":87,"y":1.515874021345053},{"x":88,"y":1.5637502040641278},{"x":89,"y":1.6095779550129374},{"x":90,"y":1.6534479800561634},{"x":91,"y":1.6953490462788887},{"x":92,"y":1.7351639432791182},{"x":93,"y":1.7726688124436725},{"x":94,"y":1.8075358783854574},{"x":95,"y":1.8393395230260274},{"x":96,"y":1.8675655478760536},{"x":97,"y":1.8916234009456208},{"x":98,"y":1.910861042246055},{"x":99,"y":1.9245820686908393},{"x":100,"y":1.9320646410224958}]}],"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50}},"value":"#gorilla_repl.vega.VegaView{:content {:axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"5ad73c2c-07da-49f8-bd98-1e4fc6917fa2\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain {:data \"5ad73c2c-07da-49f8-bd98-1e4fc6917fa2\", :field \"data.y\"}}], :marks [{:type \"line\", :from {:data \"5ad73c2c-07da-49f8-bd98-1e4fc6917fa2\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"#FF29D2\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}}], :data [{:name \"5ad73c2c-07da-49f8-bd98-1e4fc6917fa2\", :values ({:x 0, :y 2.0} {:x 1, :y 1.9284804760741014} {:x 2, :y 1.847869707654912} {:x 3, :y 1.759147940170956} {:x 4, :y 1.6633539300874196} {:x 5, :y 1.5615627444976496} {:x 6, :y 1.4548630903699578} {:x 7, :y 1.3443348754972466} {:x 8, :y 1.2310275172785075} {:x 9, :y 1.1159395524551143} {:x 10, :y 1.0000000198584484} {:x 11, :y 0.884052161868249} {:x 12, :y 0.7688397122512814} {:x 13, :y 0.6549961985048567} {:x 14, :y 0.5430374695353829} {:x 15, :y 0.4333576046967714} {:x 16, :y 0.32622828094447165} {:x 17, :y 0.2218015812237824} {:x 18, :y 0.12011612593437726} {:x 19, :y 0.02110633326264301} {:x 20, :y -0.07538547415089168} {:x 21, :y -0.16959447211977985} {:x 22, :y -0.2618166633840412} {:x 23, :y -0.35239062947183336} {:x 24, :y -0.4416781845088593} {:x 25, :y -0.5300445077444007} {:x 26, :y -0.6178383035773294} {:x 27, :y -0.705372547526059} {:x 28, :y -0.7929063568341025} {:x 29, :y -0.88062849998643} {:x 30, :y -0.9686430148760465} {:x 31, :y -1.0569573420664895} {:x 32, :y -1.1454733229161638} {:x 33, :y -1.2339813292889525} {:x 34, :y -1.3221577103472784} {:x 35, :y -1.4095656548744435} {:x 36, :y -1.4956594741653633} {:x 37, :y -1.5797922236920523} {:x 38, :y -1.6612264815843205} {:x 39, :y -1.7391480398024644} {:x 40, :y -1.8126821585831931} {:x 41, :y -1.8809119854688685} {:x 42, :y -1.9428986720936203} {:x 43, :y -1.9977026777887568} {:x 44, :y -2.044405720677371} {:x 45, :y -2.082132804456243} {:x 46, :y -2.110073758399698} {:x 47, :y -2.1275037306372786} {:x 48, :y -2.133802107707231} {:x 49, :y -2.128469353974073} {:x 50, :y -2.1111413434194244} {:x 51, :y -2.0816007929888336} {:x 52, :y -2.0397854974273417} {:x 53, :y -1.9857931381250578} {:x 54, :y -1.9198825296991564} {:x 55, :y -1.8424712574666764} {:x 56, :y -1.7541297447766149} {:x 57, :y -1.655571898771775} {:x 58, :y -1.5476425469860653} {:x 59, :y -1.431301972250347} {:x 60, :y -1.3076079261555191} {:x 61, :y -1.1776955628380277} {:x 62, :y -1.0427557895776662} {:x 63, :y -0.9040125671716012} {:x 64, :y -0.7626997233384684} {:x 65, :y -0.6200378489656382} {:x 66, :y -0.4772118474306306} {:x 67, :y -0.33534968782029867} {:x 68, :y -0.19550287980466313} {:x 69, :y -0.05862913648959707} {:x 70, :y 0.07442234839953822} {:x 71, :y 0.2029226740930839} {:x 72, :y 0.3262717872830035} {:x 73, :y 0.44400362591251114} {:x 74, :y 0.5557877774505213} {:x 75, :y 0.6614276207541374} {:x 76, :y 0.760855066199606} {:x 77, :y 0.8541220585200456} {:x 78, :y 0.9413890972582583} {:x 79, :y 1.0229111375410067} {:x 80, :y 1.0990212412755835} {:x 81, :y 1.170112478892747} {:x 82, :y 1.2366185681636452} {:x 83, :y 1.2989937911787741} {:x 84, :y 1.3576927683000244} {:x 85, :y 1.4131506108921674} {:x 86, :y 1.4657640511313044} {:x 87, :y 1.515874021345053} {:x 88, :y 1.5637502040641278} {:x 89, :y 1.6095779550129374} {:x 90, :y 1.6534479800561634} {:x 91, :y 1.6953490462788887} {:x 92, :y 1.7351639432791182} {:x 93, :y 1.7726688124436725} {:x 94, :y 1.8075358783854574} {:x 95, :y 1.8393395230260274} {:x 96, :y 1.8675655478760536} {:x 97, :y 1.8916234009456208} {:x 98, :y 1.910861042246055} {:x 99, :y 1.9245820686908393} {:x 100, :y 1.9320646410224958})}], :width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}}}"}
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
;;; 
;;; 
;;; CHOOSE THE NUMBER OF DEGREES OF FREEDOM (in this case 2 for the double prendulum)
;; **

;; @@
(def df 2)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/df</span>","value":"#'goliath/df"}
;; <=

;; @@
(def prob_inheritance 0.8) ;;to do with crossover  ;probability that terms are directly inherited from parent to child.

(def mutate_pref1 0.8)
(def mutate_pref2 0.8)
(def mutate_pref3 0.2)


; (1 - mutate_pref1) = prob of gene-delete     
; (mutate_pref1 - mutate_pref2) = prob. of gene-tweak ;;change 1 df within a gene
; (mutate_pref2 - mutate_pref3) = prob. of gene-replace,
;  mutate_pref3 = prob. of gene-add

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
      
      ;(assoc (vec (repeat length 0)) (rand-int length) 1)
    
      (new-poly-term max_pow spmax length)
      ;[0 1];;cant have this or else it will crash the mathematica
      new-gene
      )
  ))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/new-poly-term</span>","value":"#'goliath/new-poly-term"}
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
(defn gene-delete [indv] 
  (let [new-indv (subvec indv 1 (count indv))]
    (if (= new-indv []) indv
      new-indv))
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/gene-delete</span>","value":"#'goliath/gene-delete"}
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

;; **
;;; Save the population off to disk after each generation - so that possible good L can be examined later.
;;; 
;; **

;; @@
(defn SaveZ [zg str]
  (spit str (pr-str zg))
  
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/SaveZ</span>","value":"#'goliath/SaveZ"}
;; <=

;; @@
(defn task [zg]
  ;;do something
  (print (:age zg))
  (SaveZ zg "Coupled_SHO_Run1.txt")
  )
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/task</span>","value":"#'goliath/task"}
;; <=

;; @@
(def initial-zeitgeist (evolution/make-zeitgeist (random-initial-population 100 poly_length power power_sum term_length)))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/initial-zeitgeist</span>","value":"#'goliath/initial-zeitgeist"}
;; <=

;; **
;;; ## Innitial Zeitgeist Checks 
;;; 
;;; A few checks to assess the 'health' and diversity of the innitial zeitgiest population.
;;; 
;;; Considered extremely privitol in generating valid solutions, to start with a healthy population of individuals.
;;; 
;;; 
;;; **Number frequency**
;;; - The frequency the a given number appears within a genotype in a population.
;; **

;; @@
; take all poly terms from all individiuals and flatten them into one lazy sequence, perform a frequency count of the numbers present in this sequence containing all recurring numbers. 


(sort (frequencies (flatten (map :genotype (:rabble initial-zeitgeist)))))   

;Graphical plot of above.

(plot/list-plot (mapv second (sort (frequencies (flatten (map :genotype (:rabble initial-zeitgeist)))))) :colour 'red')



"Total terms of entire population :" (count (flatten (map :genotype (:rabble initial-zeitgeist))))



;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>2144</span>","value":"2144"}
;; <=

;; **
;;; **Length of Genotypes**
;;; 
;;; Lengths of generated polynomils within the innitial populaiton
;; **

;; @@
;;e.g   [[1 2 4 5] [1 2 3 4]]  length = 2

(def l_o_g (sort (frequencies (mapv count (map :genotype (:rabble initial-zeitgeist))))))

(plot/list-plot l_o_g :plot-range [:all [0  (apply max (map second l_o_g))]])

;y - frequency 
;x - length of polynomail (number of terms within a polynomial expression)
;; @@
;; =>
;;; {"type":"vega","content":{"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"1bda488e-0c8a-415a-a0d5-396448c4490d","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":[0,16]}],"marks":[{"type":"symbol","from":{"data":"1bda488e-0c8a-415a-a0d5-396448c4490d"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"fill":{"value":"steelblue"},"fillOpacity":{"value":1}},"update":{"shape":"circle","size":{"value":70},"stroke":{"value":"transparent"}},"hover":{"size":{"value":210},"stroke":{"value":"white"}}}}],"data":[{"name":"1bda488e-0c8a-415a-a0d5-396448c4490d","values":[{"x":1,"y":11},{"x":2,"y":10},{"x":3,"y":8},{"x":4,"y":9},{"x":5,"y":10},{"x":6,"y":15},{"x":7,"y":6},{"x":8,"y":16},{"x":9,"y":15}]}],"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50}},"value":"#gorilla_repl.vega.VegaView{:content {:axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"1bda488e-0c8a-415a-a0d5-396448c4490d\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain [0 16]}], :marks [{:type \"symbol\", :from {:data \"1bda488e-0c8a-415a-a0d5-396448c4490d\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :fill {:value \"steelblue\"}, :fillOpacity {:value 1}}, :update {:shape \"circle\", :size {:value 70}, :stroke {:value \"transparent\"}}, :hover {:size {:value 210}, :stroke {:value \"white\"}}}}], :data [{:name \"1bda488e-0c8a-415a-a0d5-396448c4490d\", :values ({:x 1, :y 11} {:x 2, :y 10} {:x 3, :y 8} {:x 4, :y 9} {:x 5, :y 10} {:x 6, :y 15} {:x 7, :y 6} {:x 8, :y 16} {:x 9, :y 15})}], :width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}}}"}
;; <=

;; **
;;; ** Sum of powers in a gene**
;;; 
;;; Shows the distributioons of term power sums in created polynomials. (This method does not keep track of which genes belong to which individuals.
;; **

;; @@
(defn power_sum_function
  [x]
  (apply + x))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/power_sum_function</span>","value":"#'goliath/power_sum_function"}
;; <=

;; @@
(plot/list-plot (sort (frequencies (flatten (map (fn [x] (map power_sum_function (:genotype x))) (:rabble initial-zeitgeist))))))

;y - frequency 
;x - sum of powers in polynomial term
;; @@
;; =>
;;; {"type":"vega","content":{"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"0bbfd5da-4991-4815-98b7-c346418aa324","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":{"data":"0bbfd5da-4991-4815-98b7-c346418aa324","field":"data.y"}}],"marks":[{"type":"symbol","from":{"data":"0bbfd5da-4991-4815-98b7-c346418aa324"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"fill":{"value":"steelblue"},"fillOpacity":{"value":1}},"update":{"shape":"circle","size":{"value":70},"stroke":{"value":"transparent"}},"hover":{"size":{"value":210},"stroke":{"value":"white"}}}}],"data":[{"name":"0bbfd5da-4991-4815-98b7-c346418aa324","values":[{"x":1,"y":1},{"x":2,"y":4},{"x":3,"y":15},{"x":4,"y":11},{"x":5,"y":25},{"x":6,"y":39},{"x":7,"y":68},{"x":8,"y":87},{"x":9,"y":126},{"x":10,"y":160}]}],"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50}},"value":"#gorilla_repl.vega.VegaView{:content {:axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"0bbfd5da-4991-4815-98b7-c346418aa324\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain {:data \"0bbfd5da-4991-4815-98b7-c346418aa324\", :field \"data.y\"}}], :marks [{:type \"symbol\", :from {:data \"0bbfd5da-4991-4815-98b7-c346418aa324\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :fill {:value \"steelblue\"}, :fillOpacity {:value 1}}, :update {:shape \"circle\", :size {:value 70}, :stroke {:value \"transparent\"}}, :hover {:size {:value 210}, :stroke {:value \"white\"}}}}], :data [{:name \"0bbfd5da-4991-4815-98b7-c346418aa324\", :values ({:x 1, :y 1} {:x 2, :y 4} {:x 3, :y 15} {:x 4, :y 11} {:x 5, :y 25} {:x 6, :y 39} {:x 7, :y 68} {:x 8, :y 87} {:x 9, :y 126} {:x 10, :y 160})}], :width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}}}"}
;; <=

;; @@
;;this function uses discrete differential to find speeds 
;;(goliath.mathlink.LagrangianScore/InitFunctions
  ;;mathKernelSz
  ;;"resources/"
  ;;experimentalDataSz
  ;;dt
  ;;df
  ;;)
;; @@

;; **
;;; #Reading in data
;; **

;; @@
;;this function reads in the speed accel from the initial data file
;;in the form th1, th2, th3, w1, w2, w3, a1, a2, a3
(goliath.mathlink.LagrangianScore/InitFunctionsSim
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

;; **
;;; ## Score
;; **

;; @@
(defn score
  [indv] 
	(first
      (into [] (goliath.mathlink.LagrangianScore/GetScore (into-array Integer/TYPE (vec (flatten indv))), df))
      )
  )
  
;;uncomment this to score a particular lagrangian eg coupled sho
;;(into [] (goliath.mathlink.LagrangianScore/GetScore (into-array Integer/TYPE (vec (flatten [[2 0 0 0] [0 2 0 0] [1 1 0 0] [0 0 2 0] [0 0 0 2]]))), df))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/score</span>","value":"#'goliath/score"}
;; <=

;; @@
(def memScore (memoize score))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/memScore</span>","value":"#'goliath/memScore"}
;; <=

;; **
;;; #Generation-config
;; **

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
                         :error (fn [e] (memScore e))}]
    {:ea-config              ea-config
     :score-functions        score-functions
     :reporting-function     (fn [z] (print ".") (flush))}))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/generation-config</span>","value":"#'goliath/generation-config"}
;; <=

;; @@
(time (def result (evolution/run-evolution generation-config initial-zeitgeist (fn [zg gc] (task zg) (>= (:age zg) 500)))))
;; @@
;; ->
;;; 0.1.2.3.4.5.6.7.8.9.10.11.12.13.14.15.16.17.18.19.20.21.22.23.24.25.26.27.28.29.30.31.32.33.34.35.36.37.38.39.40.41.42.43.44.45.46.47.48.49.50.51.52.53.54.55.56.57.58.59.60.61.62.63.64.65.66.67.68.69.70.71.72.73.74.75.76.77.78.79.80.81.82.83.84.85.86.87.88.89.90.91.92.93.94.95.96.97.98.99.100.101.102.103.104.105.106.107.108.109.110.111.112.113.114.115.116.117.118.119.120.121.122.123.124.125.126.
;; <-

;; @@
(goliath.mathlink.LagrangianScore/Shutdown)
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@
;((read-string (slurp "results16LScore.txt")))
;; @@

;; @@
(mapv #(println (:genotype %)) (sort-by :error (:elite result)))
;; @@
;; ->
;;; [[0 1 0 0] [1 6 0 2] [0 4 0 4] [2 1 0 4] [1 1 0 0] [1 2 1 1] [1 0 1 5] [1 6 0 0] [0 3 1 6] [1 3 1 0] [1 1 2 5]]
;;; [[0 1 0 0] [1 6 0 2] [0 4 0 4] [2 1 0 4] [1 1 0 0] [1 2 1 1] [1 0 1 5] [1 6 0 0] [0 3 1 6] [1 3 1 0] [1 1 2 5]]
;;; [[0 1 0 0] [1 6 0 2] [0 4 0 4] [2 1 0 4] [1 1 0 0] [1 2 1 1] [1 0 1 5] [1 6 0 0] [0 3 1 6] [1 3 1 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 1 0 0] [1 2 1 1] [1 0 1 5] [1 6 0 0] [0 3 1 6] [1 3 1 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 1 0 0] [1 2 1 1] [1 0 1 5] [1 6 0 0] [1 3 1 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 1 0 0] [1 2 1 1] [1 0 1 5] [1 6 0 0] [1 3 1 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 1 0 0] [1 2 1 1] [1 0 1 5] [1 6 0 0] [1 3 1 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 1 0 0] [1 2 1 1] [1 0 1 5] [1 6 0 0] [1 3 1 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 2 1 1] [1 0 1 5] [1 6 0 0] [1 3 1 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 2 1 1] [1 0 1 5] [1 6 0 0] [1 3 1 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 2 1 1] [1 0 1 5] [1 6 0 0] [1 3 1 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 2 1 1] [1 0 1 5] [1 6 0 0] [1 3 1 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 2 1 1] [1 0 1 5] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 2 1 1] [1 0 1 5] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 2 1 1] [1 0 1 5] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 2 1 1] [1 0 1 5] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 2 1 1] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 2 1 1] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 2 1 1] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 2 1 1] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [2 1 0 4] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 2] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 0]]
;;; [[0 1 0 0] [1 6 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; [[1 1 0 0]]
;;; 
;; <-
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"},{"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}],"value":"[nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil]"}
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
;;; {"type":"vega","content":{"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"f95abd38-4628-4e26-9ea7-f26b08e29bc6","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":{"data":"f95abd38-4628-4e26-9ea7-f26b08e29bc6","field":"data.y"}}],"marks":[{"type":"line","from":{"data":"f95abd38-4628-4e26-9ea7-f26b08e29bc6"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"#FF29D2"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}}],"data":[{"name":"f95abd38-4628-4e26-9ea7-f26b08e29bc6","values":[{"x":0,"y":-0.10612131636344163},{"x":1,"y":-0.12349550834839071},{"x":2,"y":-0.12349550834839071},{"x":3,"y":-0.12349550834839071},{"x":4,"y":-0.13200547798959042},{"x":5,"y":-0.13200547798959042},{"x":6,"y":-0.27599459989373276},{"x":7,"y":-0.27599459989373276},{"x":8,"y":-0.27599459989373276},{"x":9,"y":-0.27599459989373276},{"x":10,"y":-0.27599459989373276},{"x":11,"y":-0.2759945998937335},{"x":12,"y":-0.2759945998937335},{"x":13,"y":-0.2759945998937336},{"x":14,"y":-0.285083801482733},{"x":15,"y":-0.285083801482733},{"x":16,"y":-0.30840890023738304},{"x":17,"y":-0.30840890023738304},{"x":18,"y":-0.30840890023738304},{"x":19,"y":-0.30840890023738304},{"x":20,"y":-0.30840890023738304},{"x":21,"y":-0.30840890023738304},{"x":22,"y":-0.30840890023738304},{"x":23,"y":-0.30840890023738304},{"x":24,"y":-0.3085248088439997},{"x":25,"y":-0.3113212366800247},{"x":26,"y":-0.3357504690284638},{"x":27,"y":-0.3357504690284638},{"x":28,"y":-0.3357504690284638},{"x":29,"y":-0.3357504690284638},{"x":30,"y":-0.3357504690284638},{"x":31,"y":-0.3357504690284638},{"x":32,"y":-0.3357504690284638},{"x":33,"y":-0.3357504690284638},{"x":34,"y":-0.3357504690284638},{"x":35,"y":-0.33575046902846395},{"x":36,"y":-0.3357504690284641},{"x":37,"y":-0.3357504690284641},{"x":38,"y":-0.3357504690284641},{"x":39,"y":-0.3374507075136496},{"x":40,"y":-0.34311917882915244},{"x":41,"y":-0.34311917882915244},{"x":42,"y":-0.34311917882915244},{"x":43,"y":-0.34311917882915244},{"x":44,"y":-0.39067759060160295},{"x":45,"y":-0.39067759060160295},{"x":46,"y":-0.39067759060160295},{"x":47,"y":-0.3906775906016036},{"x":48,"y":-0.3906775906016036},{"x":49,"y":-0.3964185726043021},{"x":50,"y":-0.450303278041247},{"x":51,"y":-0.450303278041247},{"x":52,"y":-0.45030327804124715},{"x":53,"y":-0.45030327804124715},{"x":54,"y":-0.45030327804124715},{"x":55,"y":-0.45030327804124715},{"x":56,"y":-0.45030327804124715},{"x":57,"y":-0.45030327804124715},{"x":58,"y":-0.45030327804124715},{"x":59,"y":-0.45030327804124715},{"x":60,"y":-0.45030327804124715},{"x":61,"y":-0.4659469991466356},{"x":62,"y":-0.4659469991466356},{"x":63,"y":-0.4659469991466356},{"x":64,"y":-0.4659469991466356},{"x":65,"y":-0.4659469991466356},{"x":66,"y":-0.4659469991466356},{"x":67,"y":-0.4659469991466356},{"x":68,"y":-0.4659469991466356},{"x":69,"y":-0.4659469991466356},{"x":70,"y":-0.4659469991466356},{"x":71,"y":-0.4659469991466356},{"x":72,"y":-0.4659469991466356},{"x":73,"y":-0.4659469991466356},{"x":74,"y":-0.4784332027708211},{"x":75,"y":-0.4784332027708211},{"x":76,"y":-0.4784332027708211},{"x":77,"y":-0.4784332027708211},{"x":78,"y":-0.4784332027708211},{"x":79,"y":-0.4784332027708211},{"x":80,"y":-0.4784332027708211},{"x":81,"y":-0.48049604051417855},{"x":82,"y":-0.48049604051417855},{"x":83,"y":-0.48049604051417855},{"x":84,"y":-0.48049604051417855},{"x":85,"y":-0.48049604051417855},{"x":86,"y":-0.48049604051417855},{"x":87,"y":-0.48049604051417855},{"x":88,"y":-0.48049604051417855},{"x":89,"y":-0.48049604051417855},{"x":90,"y":-0.48049604051417855},{"x":91,"y":-0.48049604051417855},{"x":92,"y":-0.48049604051417855},{"x":93,"y":-0.48049604051417855},{"x":94,"y":-0.48049604051417855},{"x":95,"y":-0.48049604051417855},{"x":96,"y":-0.48049604051417855},{"x":97,"y":-0.48049604051417855},{"x":98,"y":-0.48049604051417855},{"x":99,"y":-0.48049604051417855},{"x":100,"y":-0.48049604051417855},{"x":101,"y":-0.48049604051417855},{"x":102,"y":-0.48049604051417855},{"x":103,"y":-0.48049604051417855},{"x":104,"y":-0.48049604051417855},{"x":105,"y":-0.48049604051417855},{"x":106,"y":-0.48049604051417855},{"x":107,"y":-0.48049604051417855},{"x":108,"y":-0.48049604051417855},{"x":109,"y":-0.48049604051417855},{"x":110,"y":-0.48049604051417855},{"x":111,"y":-0.48049604051417855},{"x":112,"y":-0.48049604051417855},{"x":113,"y":-0.5084498691339627},{"x":114,"y":-0.5084498691339627},{"x":115,"y":-0.5084498691339627},{"x":116,"y":-0.5084498691339632},{"x":117,"y":-0.5084498691339632},{"x":118,"y":-0.5084498691339632},{"x":119,"y":-0.5108212036729481},{"x":120,"y":-0.5108212036729481},{"x":121,"y":-0.5108212036729486},{"x":122,"y":-0.5108212036729486},{"x":123,"y":-0.5108212036729486},{"x":124,"y":-0.654585878059849},{"x":125,"y":-0.654585878059849},{"x":126,"y":-0.654585878059849},{"x":127,"y":-0.654585878059849},{"x":128,"y":-0.654585878059849},{"x":129,"y":-0.654585878059849},{"x":130,"y":-0.654585878059849},{"x":131,"y":-0.654585878059849},{"x":132,"y":-0.654585878059849},{"x":133,"y":-0.654585878059849},{"x":134,"y":-0.654585878059849},{"x":135,"y":-0.654585878059849},{"x":136,"y":-0.654585878059849},{"x":137,"y":-0.654585878059849},{"x":138,"y":-0.654585878059849},{"x":139,"y":-0.654585878059849},{"x":140,"y":-0.654585878059849},{"x":141,"y":-0.654585878059849},{"x":142,"y":-0.654585878059849},{"x":143,"y":-0.6589259597743009},{"x":144,"y":-0.6589259597743009},{"x":145,"y":-0.6589259597743009},{"x":146,"y":-0.6589259597743009},{"x":147,"y":-0.6589259597743009},{"x":148,"y":-0.6589259597743009},{"x":149,"y":-0.6589259597743009},{"x":150,"y":-0.6589259597743009},{"x":151,"y":-0.6589259597743009},{"x":152,"y":-0.6589259597743009},{"x":153,"y":-0.6589259597743009},{"x":154,"y":-0.6589259597743009},{"x":155,"y":-0.6589259597743009},{"x":156,"y":-0.6589259597743011},{"x":157,"y":-0.6589259597743011},{"x":158,"y":-0.6589259597743011},{"x":159,"y":-0.6589259597743011},{"x":160,"y":-0.6589259597743011},{"x":161,"y":-0.6589259597743011},{"x":162,"y":-0.6589259597743011},{"x":163,"y":-0.6589259597743011},{"x":164,"y":-0.6589259597743011},{"x":165,"y":-0.6589358064030757},{"x":166,"y":-0.6589358064030757},{"x":167,"y":-0.6589358064030765},{"x":168,"y":-0.6589358064030765},{"x":169,"y":-0.6589358064030765},{"x":170,"y":-0.6589358064030765},{"x":171,"y":-0.6632608984190425},{"x":172,"y":-0.6632608984190425},{"x":173,"y":-0.6632608984190425},{"x":174,"y":-0.6632608984190425},{"x":175,"y":-0.6632608984190425},{"x":176,"y":-0.6632608984190425},{"x":177,"y":-0.6632608984190425},{"x":178,"y":-0.6632608984190425},{"x":179,"y":-0.6632608984190425},{"x":180,"y":-0.6632608984190425},{"x":181,"y":-0.6632608984190425},{"x":182,"y":-0.6632608984190425},{"x":183,"y":-0.6632608984190425},{"x":184,"y":-0.6632608984190431},{"x":185,"y":-0.6632608984190431},{"x":186,"y":-0.6632608984190431},{"x":187,"y":-0.6632608984190431},{"x":188,"y":-0.6632608984190431},{"x":189,"y":-0.6632608984190431},{"x":190,"y":-0.6632608984190431},{"x":191,"y":-0.6632608984190431},{"x":192,"y":-0.6632608984190431},{"x":193,"y":-0.6632608984190431},{"x":194,"y":-0.6632608984190431},{"x":195,"y":-0.6632608984190431},{"x":196,"y":-0.6632608984190431},{"x":197,"y":-0.6632608984190431},{"x":198,"y":-0.6632608984190431},{"x":199,"y":-0.7232820509216807},{"x":200,"y":-0.7232820509216807},{"x":201,"y":-0.7232820509216807},{"x":202,"y":-0.7232820509216807},{"x":203,"y":-0.7232820509216807},{"x":204,"y":-0.7232820509216807},{"x":205,"y":-0.7232820509216807},{"x":206,"y":-0.7232820509216807},{"x":207,"y":-0.7232820509216807},{"x":208,"y":-0.7232820509216807},{"x":209,"y":-0.7232820509216807},{"x":210,"y":-0.7232820509216807},{"x":211,"y":-0.7232820509216807},{"x":212,"y":-0.7232820509216807},{"x":213,"y":-0.7232820509216807},{"x":214,"y":-0.7232820509216807},{"x":215,"y":-0.7232820509216807},{"x":216,"y":-0.7232820509216807},{"x":217,"y":-0.7232820509216807},{"x":218,"y":-0.7232820509216817},{"x":219,"y":-0.7232820509216817},{"x":220,"y":-0.7232820509216817},{"x":221,"y":-0.7232820509216817},{"x":222,"y":-0.7232820509216817},{"x":223,"y":-0.7237350373982939},{"x":224,"y":-0.7237350373982939},{"x":225,"y":-0.7237350373982939},{"x":226,"y":-0.7237350373982939},{"x":227,"y":-0.7237350373982949},{"x":228,"y":-0.7237350373982949},{"x":229,"y":-0.7284511912639964},{"x":230,"y":-0.7284511912639964},{"x":231,"y":-0.7284511912639964},{"x":232,"y":-0.728451191263997},{"x":233,"y":-0.728451191263997},{"x":234,"y":-0.728451191263997},{"x":235,"y":-0.728451191263997},{"x":236,"y":-0.728451191263997},{"x":237,"y":-0.728451191263997},{"x":238,"y":-0.733452807562484},{"x":239,"y":-0.733452807562484},{"x":240,"y":-0.733452807562484},{"x":241,"y":-0.7665248485724551},{"x":242,"y":-0.7665248485724551},{"x":243,"y":-0.7665248485724551},{"x":244,"y":-0.7665248485724554},{"x":245,"y":-0.7665248485724554},{"x":246,"y":-0.7665248485724554},{"x":247,"y":-0.7665248485724554},{"x":248,"y":-0.7665248485724554},{"x":249,"y":-0.7665248485724554},{"x":250,"y":-0.7665248485724554},{"x":251,"y":-0.7665248485724554},{"x":252,"y":-0.921351422602801},{"x":253,"y":-0.921351422602801},{"x":254,"y":-0.921351422602801},{"x":255,"y":-0.921351422602801},{"x":256,"y":-0.921351422602801},{"x":257,"y":-0.921351422602801},{"x":258,"y":-0.921351422602801},{"x":259,"y":-0.921351422602801},{"x":260,"y":-0.921351422602801},{"x":261,"y":-0.921351422602801},{"x":262,"y":-1.1666946303026406},{"x":263,"y":-1.1666946303026406},{"x":264,"y":-1.1666946303026406},{"x":265,"y":-1.1666946303026406},{"x":266,"y":-1.1666946303026406},{"x":267,"y":-1.1666946303026406},{"x":268,"y":-1.1666946303026406},{"x":269,"y":-1.1666946303026406},{"x":270,"y":-1.1666946303026406},{"x":271,"y":-1.1666946303026406},{"x":272,"y":-1.1666946303026406},{"x":273,"y":-1.1666946303026406},{"x":274,"y":-1.1666946303026406},{"x":275,"y":-1.1666946303026406},{"x":276,"y":-1.1666946303026406},{"x":277,"y":-1.1684211198505654},{"x":278,"y":-1.1684211198505654},{"x":279,"y":-1.1684211198505654},{"x":280,"y":-1.1684211198505654},{"x":281,"y":-1.1684211198505654},{"x":282,"y":-1.1867072412708786},{"x":283,"y":-1.3206477408608923},{"x":284,"y":-1.3206477408608923},{"x":285,"y":-1.3206477408608923},{"x":286,"y":-1.3206477408608923},{"x":287,"y":-1.3206477408608923},{"x":288,"y":-1.3206477408608956},{"x":289,"y":-1.3206477408608956},{"x":290,"y":-1.3206477408608956},{"x":291,"y":-1.3242408094773812},{"x":292,"y":-1.3242408094773812},{"x":293,"y":-1.3242408094773812},{"x":294,"y":-1.3242408094773812},{"x":295,"y":-1.3242408094773814},{"x":296,"y":-1.3242408094773814},{"x":297,"y":-1.3242408094773814},{"x":298,"y":-1.3242408094773814},{"x":299,"y":-1.3242408094773814},{"x":300,"y":-1.3242408094773828},{"x":301,"y":-1.3242408094773828},{"x":302,"y":-1.3242408094773828},{"x":303,"y":-1.3242408094773828},{"x":304,"y":-1.3323102216327334},{"x":305,"y":-1.3323102216327334},{"x":306,"y":-1.3323102216327334},{"x":307,"y":-1.3323102216327334},{"x":308,"y":-1.3323102216327334},{"x":309,"y":-1.3323102216327334},{"x":310,"y":-1.334261308714792},{"x":311,"y":-1.334261308714792},{"x":312,"y":-1.3342613087147943},{"x":313,"y":-1.3342613087147943},{"x":314,"y":-1.3342613087147943},{"x":315,"y":-1.3342613087147943},{"x":316,"y":-1.3342613087147943},{"x":317,"y":-1.3754574220302715},{"x":318,"y":-1.3754574220302715},{"x":319,"y":-1.3754574220302715},{"x":320,"y":-1.3754574220302715},{"x":321,"y":-1.3754574220302715},{"x":322,"y":-1.3754574220302715},{"x":323,"y":-1.3754574220302715},{"x":324,"y":-1.3754574220302715},{"x":325,"y":-1.3754574220302715},{"x":326,"y":-1.3754574220302715},{"x":327,"y":-1.3754574220302715},{"x":328,"y":-1.3754574220302715},{"x":329,"y":-1.3754574220302715},{"x":330,"y":-1.3754574220302715},{"x":331,"y":-1.3754574220302715},{"x":332,"y":-1.3754574220302715},{"x":333,"y":-1.3754574220302715},{"x":334,"y":-1.3754574220302715},{"x":335,"y":-1.3754574220302715},{"x":336,"y":-1.3754574220302715},{"x":337,"y":-1.3754574220302715},{"x":338,"y":-1.4347607910128237},{"x":339,"y":-1.4347607910128237},{"x":340,"y":-1.4347607910128237},{"x":341,"y":-1.4347607910128237},{"x":342,"y":-1.514169634110397},{"x":343,"y":-1.514169634110397},{"x":344,"y":-1.514169634110397},{"x":345,"y":-1.514169634110397},{"x":346,"y":-1.514169634110397},{"x":347,"y":-1.514169634110397},{"x":348,"y":-1.514169634110397},{"x":349,"y":-1.514169634110397},{"x":350,"y":-1.514169634110397},{"x":351,"y":-1.514169634110397},{"x":352,"y":-1.5167033158304333},{"x":353,"y":-1.5241034844861638},{"x":354,"y":-1.5241034844861638},{"x":355,"y":-1.5241034844861638},{"x":356,"y":-1.5241034844861638},{"x":357,"y":-1.5241034844861638},{"x":358,"y":-1.5241034844861638},{"x":359,"y":-1.5241034844861638},{"x":360,"y":-1.5241034844861638},{"x":361,"y":-1.5336072341436262},{"x":362,"y":-1.5336072341436262},{"x":363,"y":-1.5336072341436262},{"x":364,"y":-1.5336072341436262},{"x":365,"y":-1.5336072341436262},{"x":366,"y":-1.5336072341436262},{"x":367,"y":-1.5336072341436262},{"x":368,"y":-1.5336072341436262},{"x":369,"y":-1.5336072341436262},{"x":370,"y":-1.5336072341436262},{"x":371,"y":-1.5341061957679025},{"x":372,"y":-1.5341061957679025},{"x":373,"y":-1.5341061957679025},{"x":374,"y":-1.5341061957679025},{"x":375,"y":-1.5341061957679025},{"x":376,"y":-1.5341061957679025},{"x":377,"y":-1.5341061957679025},{"x":378,"y":-1.5341061957679025},{"x":379,"y":-1.5341061957679025},{"x":380,"y":-1.5341061957679025},{"x":381,"y":-1.5820654676867165},{"x":382,"y":-1.5820654676867165},{"x":383,"y":-1.5820654676867165},{"x":384,"y":-1.5820654676867165},{"x":385,"y":-1.5820654676867165},{"x":386,"y":-1.5820654676867165},{"x":387,"y":-1.5820654676867165},{"x":388,"y":-1.5820654676867165},{"x":389,"y":-1.5820654676867165},{"x":390,"y":-1.5820654676867165},{"x":391,"y":-1.5820654676867165},{"x":392,"y":-1.5820654676867165},{"x":393,"y":-1.5820654676867165},{"x":394,"y":-1.5820654676867165},{"x":395,"y":-1.5820654676867165},{"x":396,"y":-1.5820654676867165},{"x":397,"y":-1.5820654676867165},{"x":398,"y":-1.5820654676867165},{"x":399,"y":-1.5820654676867165},{"x":400,"y":-1.5820654676867165},{"x":401,"y":-1.5820654676867165},{"x":402,"y":-1.5820654676867165},{"x":403,"y":-1.5820654676867165},{"x":404,"y":-1.5911247030964386},{"x":405,"y":-1.6348756651298104},{"x":406,"y":-1.6348756651298104},{"x":407,"y":-1.6348756651298104},{"x":408,"y":-1.6348756651298104},{"x":409,"y":-1.6348756651298104},{"x":410,"y":-1.6348756651298104},{"x":411,"y":-1.6348756651298104},{"x":412,"y":-1.6348756651298104},{"x":413,"y":-1.6348756651298104},{"x":414,"y":-1.6348756651298104},{"x":415,"y":-1.6348756651298104},{"x":416,"y":-1.6348756651298104},{"x":417,"y":-1.6348756651298104},{"x":418,"y":-1.6348756651298104},{"x":419,"y":-1.6352109931652057},{"x":420,"y":-1.6352109931652057},{"x":421,"y":-1.6352109931652057},{"x":422,"y":-1.6352109931652057},{"x":423,"y":-1.636659577154861},{"x":424,"y":-1.636659577154861},{"x":425,"y":-1.636659577154861},{"x":426,"y":-1.6793586236965705},{"x":427,"y":-1.6793586236965705},{"x":428,"y":-1.6793586236965705},{"x":429,"y":-1.6793586236965705},{"x":430,"y":-1.6793586236965705},{"x":431,"y":-1.6793586236965705},{"x":432,"y":-1.6793586236965705},{"x":433,"y":-1.6793586236965723},{"x":434,"y":-1.6793586236965723},{"x":435,"y":-1.8134833955291167},{"x":436,"y":-1.925227106256439},{"x":437,"y":-1.925227106256439},{"x":438,"y":-1.925227106256439},{"x":439,"y":-1.925227106256439},{"x":440,"y":-1.925227106256439},{"x":441,"y":-1.9290925697527173},{"x":442,"y":-1.9290925697527173},{"x":443,"y":-1.9290925697527173},{"x":444,"y":-1.9290925697527173},{"x":445,"y":-1.9290925697527173},{"x":446,"y":-1.9290925697527173},{"x":447,"y":-1.9290925697527281},{"x":448,"y":-1.9421398346539256},{"x":449,"y":-1.9421398346539256},{"x":450,"y":-1.9421398346539256},{"x":451,"y":-1.9421398346539256},{"x":452,"y":-1.9421398346539256},{"x":453,"y":-1.9421398346539256},{"x":454,"y":-1.9421398346539256},{"x":455,"y":-1.9421398346539256},{"x":456,"y":-1.9721809560224968},{"x":457,"y":-1.9721809560224968},{"x":458,"y":-1.9721809560224968},{"x":459,"y":-1.9721809560224968},{"x":460,"y":-1.9721809560224968},{"x":461,"y":-1.9721809560224968},{"x":462,"y":-1.9721809560224968},{"x":463,"y":-1.9721809560225152},{"x":464,"y":-1.9721809560225152},{"x":465,"y":-1.9721809560225152},{"x":466,"y":-1.9721809560225152},{"x":467,"y":-1.9721809560225152},{"x":468,"y":-1.9721809560225152},{"x":469,"y":-1.9721809560225152},{"x":470,"y":-2.0475443733375687},{"x":471,"y":-2.0475443733375687},{"x":472,"y":-2.138704734710683},{"x":473,"y":-2.138704734710683},{"x":474,"y":-2.138704734710683},{"x":475,"y":-2.138704734710683},{"x":476,"y":-2.138704734710683},{"x":477,"y":-2.138704734710683},{"x":478,"y":-2.138704734710683},{"x":479,"y":-2.138704734710683},{"x":480,"y":-2.138704734710683},{"x":481,"y":-2.138704734710683},{"x":482,"y":-2.138704734710683},{"x":483,"y":-2.138704734710683},{"x":484,"y":-2.138704734710683},{"x":485,"y":-2.138704734710683},{"x":486,"y":-2.298761224980861},{"x":487,"y":-2.298761224980861},{"x":488,"y":-2.298761224980861},{"x":489,"y":-2.298761224980861},{"x":490,"y":-2.298761224980861},{"x":491,"y":-2.298761224980861},{"x":492,"y":-2.298761224980861},{"x":493,"y":-2.298761224980861},{"x":494,"y":-2.298761224980861},{"x":495,"y":-2.298761224980861},{"x":496,"y":-2.298761224980861},{"x":497,"y":-2.298761224980861},{"x":498,"y":-2.298761224980861},{"x":499,"y":-2.298761224980861}]}],"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50}},"value":"#gorilla_repl.vega.VegaView{:content {:axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"f95abd38-4628-4e26-9ea7-f26b08e29bc6\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain {:data \"f95abd38-4628-4e26-9ea7-f26b08e29bc6\", :field \"data.y\"}}], :marks [{:type \"line\", :from {:data \"f95abd38-4628-4e26-9ea7-f26b08e29bc6\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"#FF29D2\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}}], :data [{:name \"f95abd38-4628-4e26-9ea7-f26b08e29bc6\", :values ({:x 0, :y -0.10612131636344163} {:x 1, :y -0.12349550834839071} {:x 2, :y -0.12349550834839071} {:x 3, :y -0.12349550834839071} {:x 4, :y -0.13200547798959042} {:x 5, :y -0.13200547798959042} {:x 6, :y -0.27599459989373276} {:x 7, :y -0.27599459989373276} {:x 8, :y -0.27599459989373276} {:x 9, :y -0.27599459989373276} {:x 10, :y -0.27599459989373276} {:x 11, :y -0.2759945998937335} {:x 12, :y -0.2759945998937335} {:x 13, :y -0.2759945998937336} {:x 14, :y -0.285083801482733} {:x 15, :y -0.285083801482733} {:x 16, :y -0.30840890023738304} {:x 17, :y -0.30840890023738304} {:x 18, :y -0.30840890023738304} {:x 19, :y -0.30840890023738304} {:x 20, :y -0.30840890023738304} {:x 21, :y -0.30840890023738304} {:x 22, :y -0.30840890023738304} {:x 23, :y -0.30840890023738304} {:x 24, :y -0.3085248088439997} {:x 25, :y -0.3113212366800247} {:x 26, :y -0.3357504690284638} {:x 27, :y -0.3357504690284638} {:x 28, :y -0.3357504690284638} {:x 29, :y -0.3357504690284638} {:x 30, :y -0.3357504690284638} {:x 31, :y -0.3357504690284638} {:x 32, :y -0.3357504690284638} {:x 33, :y -0.3357504690284638} {:x 34, :y -0.3357504690284638} {:x 35, :y -0.33575046902846395} {:x 36, :y -0.3357504690284641} {:x 37, :y -0.3357504690284641} {:x 38, :y -0.3357504690284641} {:x 39, :y -0.3374507075136496} {:x 40, :y -0.34311917882915244} {:x 41, :y -0.34311917882915244} {:x 42, :y -0.34311917882915244} {:x 43, :y -0.34311917882915244} {:x 44, :y -0.39067759060160295} {:x 45, :y -0.39067759060160295} {:x 46, :y -0.39067759060160295} {:x 47, :y -0.3906775906016036} {:x 48, :y -0.3906775906016036} {:x 49, :y -0.3964185726043021} {:x 50, :y -0.450303278041247} {:x 51, :y -0.450303278041247} {:x 52, :y -0.45030327804124715} {:x 53, :y -0.45030327804124715} {:x 54, :y -0.45030327804124715} {:x 55, :y -0.45030327804124715} {:x 56, :y -0.45030327804124715} {:x 57, :y -0.45030327804124715} {:x 58, :y -0.45030327804124715} {:x 59, :y -0.45030327804124715} {:x 60, :y -0.45030327804124715} {:x 61, :y -0.4659469991466356} {:x 62, :y -0.4659469991466356} {:x 63, :y -0.4659469991466356} {:x 64, :y -0.4659469991466356} {:x 65, :y -0.4659469991466356} {:x 66, :y -0.4659469991466356} {:x 67, :y -0.4659469991466356} {:x 68, :y -0.4659469991466356} {:x 69, :y -0.4659469991466356} {:x 70, :y -0.4659469991466356} {:x 71, :y -0.4659469991466356} {:x 72, :y -0.4659469991466356} {:x 73, :y -0.4659469991466356} {:x 74, :y -0.4784332027708211} {:x 75, :y -0.4784332027708211} {:x 76, :y -0.4784332027708211} {:x 77, :y -0.4784332027708211} {:x 78, :y -0.4784332027708211} {:x 79, :y -0.4784332027708211} {:x 80, :y -0.4784332027708211} {:x 81, :y -0.48049604051417855} {:x 82, :y -0.48049604051417855} {:x 83, :y -0.48049604051417855} {:x 84, :y -0.48049604051417855} {:x 85, :y -0.48049604051417855} {:x 86, :y -0.48049604051417855} {:x 87, :y -0.48049604051417855} {:x 88, :y -0.48049604051417855} {:x 89, :y -0.48049604051417855} {:x 90, :y -0.48049604051417855} {:x 91, :y -0.48049604051417855} {:x 92, :y -0.48049604051417855} {:x 93, :y -0.48049604051417855} {:x 94, :y -0.48049604051417855} {:x 95, :y -0.48049604051417855} {:x 96, :y -0.48049604051417855} {:x 97, :y -0.48049604051417855} {:x 98, :y -0.48049604051417855} {:x 99, :y -0.48049604051417855} {:x 100, :y -0.48049604051417855} {:x 101, :y -0.48049604051417855} {:x 102, :y -0.48049604051417855} {:x 103, :y -0.48049604051417855} {:x 104, :y -0.48049604051417855} {:x 105, :y -0.48049604051417855} {:x 106, :y -0.48049604051417855} {:x 107, :y -0.48049604051417855} {:x 108, :y -0.48049604051417855} {:x 109, :y -0.48049604051417855} {:x 110, :y -0.48049604051417855} {:x 111, :y -0.48049604051417855} {:x 112, :y -0.48049604051417855} {:x 113, :y -0.5084498691339627} {:x 114, :y -0.5084498691339627} {:x 115, :y -0.5084498691339627} {:x 116, :y -0.5084498691339632} {:x 117, :y -0.5084498691339632} {:x 118, :y -0.5084498691339632} {:x 119, :y -0.5108212036729481} {:x 120, :y -0.5108212036729481} {:x 121, :y -0.5108212036729486} {:x 122, :y -0.5108212036729486} {:x 123, :y -0.5108212036729486} {:x 124, :y -0.654585878059849} {:x 125, :y -0.654585878059849} {:x 126, :y -0.654585878059849} {:x 127, :y -0.654585878059849} {:x 128, :y -0.654585878059849} {:x 129, :y -0.654585878059849} {:x 130, :y -0.654585878059849} {:x 131, :y -0.654585878059849} {:x 132, :y -0.654585878059849} {:x 133, :y -0.654585878059849} {:x 134, :y -0.654585878059849} {:x 135, :y -0.654585878059849} {:x 136, :y -0.654585878059849} {:x 137, :y -0.654585878059849} {:x 138, :y -0.654585878059849} {:x 139, :y -0.654585878059849} {:x 140, :y -0.654585878059849} {:x 141, :y -0.654585878059849} {:x 142, :y -0.654585878059849} {:x 143, :y -0.6589259597743009} {:x 144, :y -0.6589259597743009} {:x 145, :y -0.6589259597743009} {:x 146, :y -0.6589259597743009} {:x 147, :y -0.6589259597743009} {:x 148, :y -0.6589259597743009} {:x 149, :y -0.6589259597743009} {:x 150, :y -0.6589259597743009} {:x 151, :y -0.6589259597743009} {:x 152, :y -0.6589259597743009} {:x 153, :y -0.6589259597743009} {:x 154, :y -0.6589259597743009} {:x 155, :y -0.6589259597743009} {:x 156, :y -0.6589259597743011} {:x 157, :y -0.6589259597743011} {:x 158, :y -0.6589259597743011} {:x 159, :y -0.6589259597743011} {:x 160, :y -0.6589259597743011} {:x 161, :y -0.6589259597743011} {:x 162, :y -0.6589259597743011} {:x 163, :y -0.6589259597743011} {:x 164, :y -0.6589259597743011} {:x 165, :y -0.6589358064030757} {:x 166, :y -0.6589358064030757} {:x 167, :y -0.6589358064030765} {:x 168, :y -0.6589358064030765} {:x 169, :y -0.6589358064030765} {:x 170, :y -0.6589358064030765} {:x 171, :y -0.6632608984190425} {:x 172, :y -0.6632608984190425} {:x 173, :y -0.6632608984190425} {:x 174, :y -0.6632608984190425} {:x 175, :y -0.6632608984190425} {:x 176, :y -0.6632608984190425} {:x 177, :y -0.6632608984190425} {:x 178, :y -0.6632608984190425} {:x 179, :y -0.6632608984190425} {:x 180, :y -0.6632608984190425} {:x 181, :y -0.6632608984190425} {:x 182, :y -0.6632608984190425} {:x 183, :y -0.6632608984190425} {:x 184, :y -0.6632608984190431} {:x 185, :y -0.6632608984190431} {:x 186, :y -0.6632608984190431} {:x 187, :y -0.6632608984190431} {:x 188, :y -0.6632608984190431} {:x 189, :y -0.6632608984190431} {:x 190, :y -0.6632608984190431} {:x 191, :y -0.6632608984190431} {:x 192, :y -0.6632608984190431} {:x 193, :y -0.6632608984190431} {:x 194, :y -0.6632608984190431} {:x 195, :y -0.6632608984190431} {:x 196, :y -0.6632608984190431} {:x 197, :y -0.6632608984190431} {:x 198, :y -0.6632608984190431} {:x 199, :y -0.7232820509216807} {:x 200, :y -0.7232820509216807} {:x 201, :y -0.7232820509216807} {:x 202, :y -0.7232820509216807} {:x 203, :y -0.7232820509216807} {:x 204, :y -0.7232820509216807} {:x 205, :y -0.7232820509216807} {:x 206, :y -0.7232820509216807} {:x 207, :y -0.7232820509216807} {:x 208, :y -0.7232820509216807} {:x 209, :y -0.7232820509216807} {:x 210, :y -0.7232820509216807} {:x 211, :y -0.7232820509216807} {:x 212, :y -0.7232820509216807} {:x 213, :y -0.7232820509216807} {:x 214, :y -0.7232820509216807} {:x 215, :y -0.7232820509216807} {:x 216, :y -0.7232820509216807} {:x 217, :y -0.7232820509216807} {:x 218, :y -0.7232820509216817} {:x 219, :y -0.7232820509216817} {:x 220, :y -0.7232820509216817} {:x 221, :y -0.7232820509216817} {:x 222, :y -0.7232820509216817} {:x 223, :y -0.7237350373982939} {:x 224, :y -0.7237350373982939} {:x 225, :y -0.7237350373982939} {:x 226, :y -0.7237350373982939} {:x 227, :y -0.7237350373982949} {:x 228, :y -0.7237350373982949} {:x 229, :y -0.7284511912639964} {:x 230, :y -0.7284511912639964} {:x 231, :y -0.7284511912639964} {:x 232, :y -0.728451191263997} {:x 233, :y -0.728451191263997} {:x 234, :y -0.728451191263997} {:x 235, :y -0.728451191263997} {:x 236, :y -0.728451191263997} {:x 237, :y -0.728451191263997} {:x 238, :y -0.733452807562484} {:x 239, :y -0.733452807562484} {:x 240, :y -0.733452807562484} {:x 241, :y -0.7665248485724551} {:x 242, :y -0.7665248485724551} {:x 243, :y -0.7665248485724551} {:x 244, :y -0.7665248485724554} {:x 245, :y -0.7665248485724554} {:x 246, :y -0.7665248485724554} {:x 247, :y -0.7665248485724554} {:x 248, :y -0.7665248485724554} {:x 249, :y -0.7665248485724554} {:x 250, :y -0.7665248485724554} {:x 251, :y -0.7665248485724554} {:x 252, :y -0.921351422602801} {:x 253, :y -0.921351422602801} {:x 254, :y -0.921351422602801} {:x 255, :y -0.921351422602801} {:x 256, :y -0.921351422602801} {:x 257, :y -0.921351422602801} {:x 258, :y -0.921351422602801} {:x 259, :y -0.921351422602801} {:x 260, :y -0.921351422602801} {:x 261, :y -0.921351422602801} {:x 262, :y -1.1666946303026406} {:x 263, :y -1.1666946303026406} {:x 264, :y -1.1666946303026406} {:x 265, :y -1.1666946303026406} {:x 266, :y -1.1666946303026406} {:x 267, :y -1.1666946303026406} {:x 268, :y -1.1666946303026406} {:x 269, :y -1.1666946303026406} {:x 270, :y -1.1666946303026406} {:x 271, :y -1.1666946303026406} {:x 272, :y -1.1666946303026406} {:x 273, :y -1.1666946303026406} {:x 274, :y -1.1666946303026406} {:x 275, :y -1.1666946303026406} {:x 276, :y -1.1666946303026406} {:x 277, :y -1.1684211198505654} {:x 278, :y -1.1684211198505654} {:x 279, :y -1.1684211198505654} {:x 280, :y -1.1684211198505654} {:x 281, :y -1.1684211198505654} {:x 282, :y -1.1867072412708786} {:x 283, :y -1.3206477408608923} {:x 284, :y -1.3206477408608923} {:x 285, :y -1.3206477408608923} {:x 286, :y -1.3206477408608923} {:x 287, :y -1.3206477408608923} {:x 288, :y -1.3206477408608956} {:x 289, :y -1.3206477408608956} {:x 290, :y -1.3206477408608956} {:x 291, :y -1.3242408094773812} {:x 292, :y -1.3242408094773812} {:x 293, :y -1.3242408094773812} {:x 294, :y -1.3242408094773812} {:x 295, :y -1.3242408094773814} {:x 296, :y -1.3242408094773814} {:x 297, :y -1.3242408094773814} {:x 298, :y -1.3242408094773814} {:x 299, :y -1.3242408094773814} {:x 300, :y -1.3242408094773828} {:x 301, :y -1.3242408094773828} {:x 302, :y -1.3242408094773828} {:x 303, :y -1.3242408094773828} {:x 304, :y -1.3323102216327334} {:x 305, :y -1.3323102216327334} {:x 306, :y -1.3323102216327334} {:x 307, :y -1.3323102216327334} {:x 308, :y -1.3323102216327334} {:x 309, :y -1.3323102216327334} {:x 310, :y -1.334261308714792} {:x 311, :y -1.334261308714792} {:x 312, :y -1.3342613087147943} {:x 313, :y -1.3342613087147943} {:x 314, :y -1.3342613087147943} {:x 315, :y -1.3342613087147943} {:x 316, :y -1.3342613087147943} {:x 317, :y -1.3754574220302715} {:x 318, :y -1.3754574220302715} {:x 319, :y -1.3754574220302715} {:x 320, :y -1.3754574220302715} {:x 321, :y -1.3754574220302715} {:x 322, :y -1.3754574220302715} {:x 323, :y -1.3754574220302715} {:x 324, :y -1.3754574220302715} {:x 325, :y -1.3754574220302715} {:x 326, :y -1.3754574220302715} {:x 327, :y -1.3754574220302715} {:x 328, :y -1.3754574220302715} {:x 329, :y -1.3754574220302715} {:x 330, :y -1.3754574220302715} {:x 331, :y -1.3754574220302715} {:x 332, :y -1.3754574220302715} {:x 333, :y -1.3754574220302715} {:x 334, :y -1.3754574220302715} {:x 335, :y -1.3754574220302715} {:x 336, :y -1.3754574220302715} {:x 337, :y -1.3754574220302715} {:x 338, :y -1.4347607910128237} {:x 339, :y -1.4347607910128237} {:x 340, :y -1.4347607910128237} {:x 341, :y -1.4347607910128237} {:x 342, :y -1.514169634110397} {:x 343, :y -1.514169634110397} {:x 344, :y -1.514169634110397} {:x 345, :y -1.514169634110397} {:x 346, :y -1.514169634110397} {:x 347, :y -1.514169634110397} {:x 348, :y -1.514169634110397} {:x 349, :y -1.514169634110397} {:x 350, :y -1.514169634110397} {:x 351, :y -1.514169634110397} {:x 352, :y -1.5167033158304333} {:x 353, :y -1.5241034844861638} {:x 354, :y -1.5241034844861638} {:x 355, :y -1.5241034844861638} {:x 356, :y -1.5241034844861638} {:x 357, :y -1.5241034844861638} {:x 358, :y -1.5241034844861638} {:x 359, :y -1.5241034844861638} {:x 360, :y -1.5241034844861638} {:x 361, :y -1.5336072341436262} {:x 362, :y -1.5336072341436262} {:x 363, :y -1.5336072341436262} {:x 364, :y -1.5336072341436262} {:x 365, :y -1.5336072341436262} {:x 366, :y -1.5336072341436262} {:x 367, :y -1.5336072341436262} {:x 368, :y -1.5336072341436262} {:x 369, :y -1.5336072341436262} {:x 370, :y -1.5336072341436262} {:x 371, :y -1.5341061957679025} {:x 372, :y -1.5341061957679025} {:x 373, :y -1.5341061957679025} {:x 374, :y -1.5341061957679025} {:x 375, :y -1.5341061957679025} {:x 376, :y -1.5341061957679025} {:x 377, :y -1.5341061957679025} {:x 378, :y -1.5341061957679025} {:x 379, :y -1.5341061957679025} {:x 380, :y -1.5341061957679025} {:x 381, :y -1.5820654676867165} {:x 382, :y -1.5820654676867165} {:x 383, :y -1.5820654676867165} {:x 384, :y -1.5820654676867165} {:x 385, :y -1.5820654676867165} {:x 386, :y -1.5820654676867165} {:x 387, :y -1.5820654676867165} {:x 388, :y -1.5820654676867165} {:x 389, :y -1.5820654676867165} {:x 390, :y -1.5820654676867165} {:x 391, :y -1.5820654676867165} {:x 392, :y -1.5820654676867165} {:x 393, :y -1.5820654676867165} {:x 394, :y -1.5820654676867165} {:x 395, :y -1.5820654676867165} {:x 396, :y -1.5820654676867165} {:x 397, :y -1.5820654676867165} {:x 398, :y -1.5820654676867165} {:x 399, :y -1.5820654676867165} {:x 400, :y -1.5820654676867165} {:x 401, :y -1.5820654676867165} {:x 402, :y -1.5820654676867165} {:x 403, :y -1.5820654676867165} {:x 404, :y -1.5911247030964386} {:x 405, :y -1.6348756651298104} {:x 406, :y -1.6348756651298104} {:x 407, :y -1.6348756651298104} {:x 408, :y -1.6348756651298104} {:x 409, :y -1.6348756651298104} {:x 410, :y -1.6348756651298104} {:x 411, :y -1.6348756651298104} {:x 412, :y -1.6348756651298104} {:x 413, :y -1.6348756651298104} {:x 414, :y -1.6348756651298104} {:x 415, :y -1.6348756651298104} {:x 416, :y -1.6348756651298104} {:x 417, :y -1.6348756651298104} {:x 418, :y -1.6348756651298104} {:x 419, :y -1.6352109931652057} {:x 420, :y -1.6352109931652057} {:x 421, :y -1.6352109931652057} {:x 422, :y -1.6352109931652057} {:x 423, :y -1.636659577154861} {:x 424, :y -1.636659577154861} {:x 425, :y -1.636659577154861} {:x 426, :y -1.6793586236965705} {:x 427, :y -1.6793586236965705} {:x 428, :y -1.6793586236965705} {:x 429, :y -1.6793586236965705} {:x 430, :y -1.6793586236965705} {:x 431, :y -1.6793586236965705} {:x 432, :y -1.6793586236965705} {:x 433, :y -1.6793586236965723} {:x 434, :y -1.6793586236965723} {:x 435, :y -1.8134833955291167} {:x 436, :y -1.925227106256439} {:x 437, :y -1.925227106256439} {:x 438, :y -1.925227106256439} {:x 439, :y -1.925227106256439} {:x 440, :y -1.925227106256439} {:x 441, :y -1.9290925697527173} {:x 442, :y -1.9290925697527173} {:x 443, :y -1.9290925697527173} {:x 444, :y -1.9290925697527173} {:x 445, :y -1.9290925697527173} {:x 446, :y -1.9290925697527173} {:x 447, :y -1.9290925697527281} {:x 448, :y -1.9421398346539256} {:x 449, :y -1.9421398346539256} {:x 450, :y -1.9421398346539256} {:x 451, :y -1.9421398346539256} {:x 452, :y -1.9421398346539256} {:x 453, :y -1.9421398346539256} {:x 454, :y -1.9421398346539256} {:x 455, :y -1.9421398346539256} {:x 456, :y -1.9721809560224968} {:x 457, :y -1.9721809560224968} {:x 458, :y -1.9721809560224968} {:x 459, :y -1.9721809560224968} {:x 460, :y -1.9721809560224968} {:x 461, :y -1.9721809560224968} {:x 462, :y -1.9721809560224968} {:x 463, :y -1.9721809560225152} {:x 464, :y -1.9721809560225152} {:x 465, :y -1.9721809560225152} {:x 466, :y -1.9721809560225152} {:x 467, :y -1.9721809560225152} {:x 468, :y -1.9721809560225152} {:x 469, :y -1.9721809560225152} {:x 470, :y -2.0475443733375687} {:x 471, :y -2.0475443733375687} {:x 472, :y -2.138704734710683} {:x 473, :y -2.138704734710683} {:x 474, :y -2.138704734710683} {:x 475, :y -2.138704734710683} {:x 476, :y -2.138704734710683} {:x 477, :y -2.138704734710683} {:x 478, :y -2.138704734710683} {:x 479, :y -2.138704734710683} {:x 480, :y -2.138704734710683} {:x 481, :y -2.138704734710683} {:x 482, :y -2.138704734710683} {:x 483, :y -2.138704734710683} {:x 484, :y -2.138704734710683} {:x 485, :y -2.138704734710683} {:x 486, :y -2.298761224980861} {:x 487, :y -2.298761224980861} {:x 488, :y -2.298761224980861} {:x 489, :y -2.298761224980861} {:x 490, :y -2.298761224980861} {:x 491, :y -2.298761224980861} {:x 492, :y -2.298761224980861} {:x 493, :y -2.298761224980861} {:x 494, :y -2.298761224980861} {:x 495, :y -2.298761224980861} {:x 496, :y -2.298761224980861} {:x 497, :y -2.298761224980861} {:x 498, :y -2.298761224980861} {:x 499, :y -2.298761224980861})}], :width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}}}"}
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
;;; {"type":"vega","content":{"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50},"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"85f99bf7-f834-4b0f-bfab-4578e9502c38","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":{"data":"85f99bf7-f834-4b0f-bfab-4578e9502c38","field":"data.y"}}],"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"data":[{"name":"85f99bf7-f834-4b0f-bfab-4578e9502c38","values":[{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-0.8059581442333353,"y":4},{"x":-0.4589393529282559,"y":2},{"x":1.4368469964806807,"y":1},{"x":-0.4589393529282559,"y":2},{"x":1.4368469964806807,"y":1},{"x":-0.4589393529282559,"y":2},{"x":-0.4589393529282559,"y":2},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-0.5942926545724515,"y":3},{"x":-0.4589393529282559,"y":2},{"x":-0.8059581442333353,"y":4},{"x":-0.8059581442333353,"y":4},{"x":-0.4589393529282559,"y":2},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-0.4589393529282559,"y":2},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-0.5942926545724515,"y":3},{"x":1.4368469964806807,"y":1},{"x":-0.5942926545724515,"y":3},{"x":1.4368469964806807,"y":1},{"x":-0.4589393529282559,"y":2},{"x":1.4368469964806807,"y":1},{"x":-0.4589393529282559,"y":2},{"x":-0.4589393529282559,"y":2},{"x":-0.5942926545724515,"y":3},{"x":1.4368469964806807,"y":1},{"x":-0.4589393529282559,"y":2},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-0.5942926545724515,"y":3},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-0.5942926545724515,"y":3},{"x":-0.8059581442333353,"y":4},{"x":1.4368469964806807,"y":1},{"x":-0.8059581442333353,"y":4},{"x":-0.4589393529282559,"y":2},{"x":-0.4589393529282559,"y":2},{"x":-0.5942926545724515,"y":3},{"x":1.4368469964806807,"y":1},{"x":-0.4589393529282559,"y":2},{"x":1.4368469964806807,"y":1},{"x":-0.8059581442333353,"y":4},{"x":-0.5942926545724515,"y":3},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-0.4589393529282559,"y":2},{"x":-0.5942926545724515,"y":3},{"x":1.4368469964806807,"y":1},{"x":-0.4589393529282559,"y":2},{"x":-0.8059581442333353,"y":4},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-0.4589393529282559,"y":2},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-0.4589393529282559,"y":2},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-0.5942926545724515,"y":3},{"x":1.4368469964806807,"y":1},{"x":-0.4589393529282559,"y":2},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-0.8059581442333353,"y":4},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-2.0072115684121643,"y":8},{"x":-2.0072115684121643,"y":8},{"x":-2.0072115684121643,"y":8},{"x":-2.0072115684121643,"y":8},{"x":-1.9099417631646933,"y":7},{"x":-1.9099417631646933,"y":7},{"x":-1.9099417631646933,"y":7},{"x":-1.9099417631646933,"y":7},{"x":-1.7374705404343287,"y":6},{"x":-1.7374705404343287,"y":6},{"x":-1.7374705404343287,"y":6},{"x":-1.7374705404343287,"y":6},{"x":-1.452985831073355,"y":5},{"x":-1.452985831073355,"y":5},{"x":-1.452985831073355,"y":5},{"x":-1.452985831073355,"y":5},{"x":-2.298761224980861,"y":11},{"x":-2.298761224980861,"y":11},{"x":-2.17727991716631,"y":10},{"x":-2.130653782248272,"y":9}]},{"name":"727d244f-8a9a-4fed-b01d-7404710b05ce","values":[{"x":7.447802463243995,"y":4},{"x":16.212236546842988,"y":1},{"x":-1.9124537837571551,"y":8},{"x":-2.0485529509367697,"y":9},{"x":-0.8499541930829045,"y":6},{"x":7.447802463243995,"y":4},{"x":-0.16513265665558574,"y":8},{"x":16.212236546842988,"y":1},{"x":-0.6725656968989567,"y":5},{"x":-0.8708924736711212,"y":10},{"x":-0.047325776291742075,"y":4},{"x":-0.4589393529282559,"y":2},{"x":-0.5878150141385403,"y":5},{"x":6.403828512994876,"y":5},{"x":-0.047325776291742075,"y":4},{"x":1.37311423617054,"y":10},{"x":1.4368469964806807,"y":1},{"x":-0.14434015375772213,"y":5},{"x":-0.6725656968989567,"y":5},{"x":-0.4589393529282559,"y":2},{"x":8.39513969183821,"y":1},{"x":-0.5942926545724515,"y":3},{"x":-0.5529145992705946,"y":4},{"x":7.723912758363597,"y":2},{"x":1.4260272213542324,"y":4},{"x":-0.8380237774005208,"y":6},{"x":-0.047325776291741846,"y":3},{"x":8.18490226752753,"y":2},{"x":-0.5529145992705946,"y":4},{"x":-0.6949785389580226,"y":5},{"x":-0.047325776291741846,"y":3},{"x":1.4298569386607616,"y":2},{"x":1.4298569386607616,"y":2},{"x":-0.04634141204915413,"y":2},{"x":-0.559895414085866,"y":5},{"x":-0.5942926545724515,"y":3},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-0.5942926545724515,"y":3},{"x":-0.4589393529282559,"y":2},{"x":-0.5942926545724515,"y":3},{"x":1.4368469964806807,"y":1},{"x":-1.7374705404343287,"y":6},{"x":1.4368469964806807,"y":1},{"x":-0.5881461352194489,"y":7},{"x":1.4012469977291386,"y":4},{"x":-0.8059581442333353,"y":4},{"x":-0.5942926545724515,"y":3},{"x":-1.4941471160678161,"y":6},{"x":-0.1336509489717529,"y":6},{"x":-0.5942926545724515,"y":3},{"x":-0.8059581442333353,"y":4},{"x":1.419428813816475,"y":2},{"x":-0.5471440601847106,"y":3},{"x":1.42977250522492,"y":3},{"x":-0.04634141204915413,"y":2},{"x":-0.13509422429276782,"y":6},{"x":-1.7374705404343287,"y":6},{"x":-0.8455397130194695,"y":5},{"x":-1.9099417631646933,"y":7},{"x":-1.7374705404343287,"y":6},{"x":-0.8982685559359308,"y":6},{"x":-0.4589393529282559,"y":2},{"x":-0.047325776291741846,"y":3},{"x":6.9733647602425926,"y":5},{"x":-0.4891436964557815,"y":3},{"x":7.447802463243995,"y":3},{"x":-0.6949785389580226,"y":5},{"x":-0.9555286355244197,"y":7},{"x":-0.5529145992705946,"y":4},{"x":7.723912758363597,"y":2},{"x":-0.0878688352963173,"y":2},{"x":-0.46155563912247133,"y":3},{"x":-0.13362794882387027,"y":4},{"x":-1.452985831073355,"y":5},{"x":-0.04634141204915413,"y":2},{"x":-0.15197678660757977,"y":3},{"x":7.447802463243995,"y":3},{"x":16.378481015015257,"y":1},{"x":-0.4891436964557815,"y":3},{"x":-0.6083963424280937,"y":4},{"x":-1.452985831073355,"y":5},{"x":-0.8499541930829051,"y":5},{"x":1.4122894282888883,"y":3},{"x":-0.5942926545724515,"y":3},{"x":1.4368469964806807,"y":1},{"x":-0.5364118590365204,"y":5},{"x":6.403828512994876,"y":4},{"x":-1.7374705404343287,"y":6},{"x":-1.452985831073355,"y":5},{"x":1.436152807470486,"y":2},{"x":-0.9555286355244197,"y":7},{"x":-0.4589393529282559,"y":2},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-2.130653782248272,"y":9},{"x":16.378481015015257,"y":1},{"x":-2.0072115684121643,"y":8},{"x":-0.5529145992705946,"y":4},{"x":-0.04667276545130464,"y":4}]},{"name":"e3bcad74-936b-4b40-a31a-5fdf5d03a8c1","values":[{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-0.8059581442333353,"y":4},{"x":-0.4589393529282559,"y":2},{"x":1.4368469964806807,"y":1},{"x":-0.4589393529282559,"y":2},{"x":1.4368469964806807,"y":1},{"x":-0.4589393529282559,"y":2},{"x":-0.4589393529282559,"y":2},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-0.5942926545724515,"y":3},{"x":-0.4589393529282559,"y":2},{"x":-0.8059581442333353,"y":4},{"x":-0.8059581442333353,"y":4},{"x":-0.4589393529282559,"y":2},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-0.4589393529282559,"y":2},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-0.5942926545724515,"y":3},{"x":1.4368469964806807,"y":1},{"x":-0.5942926545724515,"y":3},{"x":1.4368469964806807,"y":1},{"x":-0.4589393529282559,"y":2},{"x":1.4368469964806807,"y":1},{"x":-0.4589393529282559,"y":2},{"x":-0.4589393529282559,"y":2},{"x":-0.5942926545724515,"y":3},{"x":1.4368469964806807,"y":1},{"x":-0.4589393529282559,"y":2},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-0.5942926545724515,"y":3},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-0.5942926545724515,"y":3},{"x":-0.8059581442333353,"y":4},{"x":1.4368469964806807,"y":1},{"x":-0.8059581442333353,"y":4},{"x":-0.4589393529282559,"y":2},{"x":-0.4589393529282559,"y":2},{"x":-0.5942926545724515,"y":3},{"x":1.4368469964806807,"y":1},{"x":-0.4589393529282559,"y":2},{"x":1.4368469964806807,"y":1},{"x":-0.8059581442333353,"y":4},{"x":-0.5942926545724515,"y":3},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-0.4589393529282559,"y":2},{"x":-0.5942926545724515,"y":3},{"x":1.4368469964806807,"y":1},{"x":-0.4589393529282559,"y":2},{"x":-0.8059581442333353,"y":4},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-0.4589393529282559,"y":2},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-0.4589393529282559,"y":2},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-0.5942926545724515,"y":3},{"x":1.4368469964806807,"y":1},{"x":-0.4589393529282559,"y":2},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-0.8059581442333353,"y":4},{"x":1.4368469964806807,"y":1},{"x":1.4368469964806807,"y":1},{"x":-2.0072115684121643,"y":8},{"x":-2.0072115684121643,"y":8},{"x":-2.0072115684121643,"y":8},{"x":-2.0072115684121643,"y":8},{"x":-1.9099417631646933,"y":7},{"x":-1.9099417631646933,"y":7},{"x":-1.9099417631646933,"y":7},{"x":-1.9099417631646933,"y":7},{"x":-1.7374705404343287,"y":6},{"x":-1.7374705404343287,"y":6},{"x":-1.7374705404343287,"y":6},{"x":-1.7374705404343287,"y":6},{"x":-1.452985831073355,"y":5},{"x":-1.452985831073355,"y":5},{"x":-1.452985831073355,"y":5},{"x":-1.452985831073355,"y":5},{"x":-2.298761224980861,"y":11},{"x":-2.298761224980861,"y":11},{"x":-2.17727991716631,"y":10},{"x":-2.130653782248272,"y":9}]}],"marks":[{"type":"symbol","from":{"data":"85f99bf7-f834-4b0f-bfab-4578e9502c38"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"fill":{"value":"red"},"fillOpacity":{"value":1}},"update":{"shape":"circle","size":{"value":70},"stroke":{"value":"transparent"}},"hover":{"size":{"value":210},"stroke":{"value":"white"}}}},{"type":"symbol","from":{"data":"727d244f-8a9a-4fed-b01d-7404710b05ce"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"fill":{"value":"blue"},"fillOpacity":{"value":1}},"update":{"shape":"circle","size":{"value":70},"stroke":{"value":"transparent"}},"hover":{"size":{"value":210},"stroke":{"value":"white"}}}},{"type":"symbol","from":{"data":"e3bcad74-936b-4b40-a31a-5fdf5d03a8c1"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"fill":{"value":"#ff29d2"},"fillOpacity":{"value":1}},"update":{"shape":"circle","size":{"value":70},"stroke":{"value":"transparent"}},"hover":{"size":{"value":210},"stroke":{"value":"white"}}}}]},"value":"#gorilla_repl.vega.VegaView{:content {:width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}, :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"85f99bf7-f834-4b0f-bfab-4578e9502c38\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain {:data \"85f99bf7-f834-4b0f-bfab-4578e9502c38\", :field \"data.y\"}}], :axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :data ({:name \"85f99bf7-f834-4b0f-bfab-4578e9502c38\", :values ({:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -0.8059581442333353, :y 4} {:x -0.4589393529282559, :y 2} {:x 1.4368469964806807, :y 1} {:x -0.4589393529282559, :y 2} {:x 1.4368469964806807, :y 1} {:x -0.4589393529282559, :y 2} {:x -0.4589393529282559, :y 2} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -0.5942926545724515, :y 3} {:x -0.4589393529282559, :y 2} {:x -0.8059581442333353, :y 4} {:x -0.8059581442333353, :y 4} {:x -0.4589393529282559, :y 2} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -0.4589393529282559, :y 2} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -0.5942926545724515, :y 3} {:x 1.4368469964806807, :y 1} {:x -0.5942926545724515, :y 3} {:x 1.4368469964806807, :y 1} {:x -0.4589393529282559, :y 2} {:x 1.4368469964806807, :y 1} {:x -0.4589393529282559, :y 2} {:x -0.4589393529282559, :y 2} {:x -0.5942926545724515, :y 3} {:x 1.4368469964806807, :y 1} {:x -0.4589393529282559, :y 2} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -0.5942926545724515, :y 3} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -0.5942926545724515, :y 3} {:x -0.8059581442333353, :y 4} {:x 1.4368469964806807, :y 1} {:x -0.8059581442333353, :y 4} {:x -0.4589393529282559, :y 2} {:x -0.4589393529282559, :y 2} {:x -0.5942926545724515, :y 3} {:x 1.4368469964806807, :y 1} {:x -0.4589393529282559, :y 2} {:x 1.4368469964806807, :y 1} {:x -0.8059581442333353, :y 4} {:x -0.5942926545724515, :y 3} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -0.4589393529282559, :y 2} {:x -0.5942926545724515, :y 3} {:x 1.4368469964806807, :y 1} {:x -0.4589393529282559, :y 2} {:x -0.8059581442333353, :y 4} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -0.4589393529282559, :y 2} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -0.4589393529282559, :y 2} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -0.5942926545724515, :y 3} {:x 1.4368469964806807, :y 1} {:x -0.4589393529282559, :y 2} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -0.8059581442333353, :y 4} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -2.0072115684121643, :y 8} {:x -2.0072115684121643, :y 8} {:x -2.0072115684121643, :y 8} {:x -2.0072115684121643, :y 8} {:x -1.9099417631646933, :y 7} {:x -1.9099417631646933, :y 7} {:x -1.9099417631646933, :y 7} {:x -1.9099417631646933, :y 7} {:x -1.7374705404343287, :y 6} {:x -1.7374705404343287, :y 6} {:x -1.7374705404343287, :y 6} {:x -1.7374705404343287, :y 6} {:x -1.452985831073355, :y 5} {:x -1.452985831073355, :y 5} {:x -1.452985831073355, :y 5} {:x -1.452985831073355, :y 5} {:x -2.298761224980861, :y 11} {:x -2.298761224980861, :y 11} {:x -2.17727991716631, :y 10} {:x -2.130653782248272, :y 9})} {:name \"727d244f-8a9a-4fed-b01d-7404710b05ce\", :values ({:x 7.447802463243995, :y 4} {:x 16.212236546842988, :y 1} {:x -1.9124537837571551, :y 8} {:x -2.0485529509367697, :y 9} {:x -0.8499541930829045, :y 6} {:x 7.447802463243995, :y 4} {:x -0.16513265665558574, :y 8} {:x 16.212236546842988, :y 1} {:x -0.6725656968989567, :y 5} {:x -0.8708924736711212, :y 10} {:x -0.047325776291742075, :y 4} {:x -0.4589393529282559, :y 2} {:x -0.5878150141385403, :y 5} {:x 6.403828512994876, :y 5} {:x -0.047325776291742075, :y 4} {:x 1.37311423617054, :y 10} {:x 1.4368469964806807, :y 1} {:x -0.14434015375772213, :y 5} {:x -0.6725656968989567, :y 5} {:x -0.4589393529282559, :y 2} {:x 8.39513969183821, :y 1} {:x -0.5942926545724515, :y 3} {:x -0.5529145992705946, :y 4} {:x 7.723912758363597, :y 2} {:x 1.4260272213542324, :y 4} {:x -0.8380237774005208, :y 6} {:x -0.047325776291741846, :y 3} {:x 8.18490226752753, :y 2} {:x -0.5529145992705946, :y 4} {:x -0.6949785389580226, :y 5} {:x -0.047325776291741846, :y 3} {:x 1.4298569386607616, :y 2} {:x 1.4298569386607616, :y 2} {:x -0.04634141204915413, :y 2} {:x -0.559895414085866, :y 5} {:x -0.5942926545724515, :y 3} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -0.5942926545724515, :y 3} {:x -0.4589393529282559, :y 2} {:x -0.5942926545724515, :y 3} {:x 1.4368469964806807, :y 1} {:x -1.7374705404343287, :y 6} {:x 1.4368469964806807, :y 1} {:x -0.5881461352194489, :y 7} {:x 1.4012469977291386, :y 4} {:x -0.8059581442333353, :y 4} {:x -0.5942926545724515, :y 3} {:x -1.4941471160678161, :y 6} {:x -0.1336509489717529, :y 6} {:x -0.5942926545724515, :y 3} {:x -0.8059581442333353, :y 4} {:x 1.419428813816475, :y 2} {:x -0.5471440601847106, :y 3} {:x 1.42977250522492, :y 3} {:x -0.04634141204915413, :y 2} {:x -0.13509422429276782, :y 6} {:x -1.7374705404343287, :y 6} {:x -0.8455397130194695, :y 5} {:x -1.9099417631646933, :y 7} {:x -1.7374705404343287, :y 6} {:x -0.8982685559359308, :y 6} {:x -0.4589393529282559, :y 2} {:x -0.047325776291741846, :y 3} {:x 6.9733647602425926, :y 5} {:x -0.4891436964557815, :y 3} {:x 7.447802463243995, :y 3} {:x -0.6949785389580226, :y 5} {:x -0.9555286355244197, :y 7} {:x -0.5529145992705946, :y 4} {:x 7.723912758363597, :y 2} {:x -0.0878688352963173, :y 2} {:x -0.46155563912247133, :y 3} {:x -0.13362794882387027, :y 4} {:x -1.452985831073355, :y 5} {:x -0.04634141204915413, :y 2} {:x -0.15197678660757977, :y 3} {:x 7.447802463243995, :y 3} {:x 16.378481015015257, :y 1} {:x -0.4891436964557815, :y 3} {:x -0.6083963424280937, :y 4} {:x -1.452985831073355, :y 5} {:x -0.8499541930829051, :y 5} {:x 1.4122894282888883, :y 3} {:x -0.5942926545724515, :y 3} {:x 1.4368469964806807, :y 1} {:x -0.5364118590365204, :y 5} {:x 6.403828512994876, :y 4} {:x -1.7374705404343287, :y 6} {:x -1.452985831073355, :y 5} {:x 1.436152807470486, :y 2} {:x -0.9555286355244197, :y 7} {:x -0.4589393529282559, :y 2} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -2.130653782248272, :y 9} {:x 16.378481015015257, :y 1} {:x -2.0072115684121643, :y 8} {:x -0.5529145992705946, :y 4} {:x -0.04667276545130464, :y 4})} {:name \"e3bcad74-936b-4b40-a31a-5fdf5d03a8c1\", :values ({:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -0.8059581442333353, :y 4} {:x -0.4589393529282559, :y 2} {:x 1.4368469964806807, :y 1} {:x -0.4589393529282559, :y 2} {:x 1.4368469964806807, :y 1} {:x -0.4589393529282559, :y 2} {:x -0.4589393529282559, :y 2} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -0.5942926545724515, :y 3} {:x -0.4589393529282559, :y 2} {:x -0.8059581442333353, :y 4} {:x -0.8059581442333353, :y 4} {:x -0.4589393529282559, :y 2} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -0.4589393529282559, :y 2} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -0.5942926545724515, :y 3} {:x 1.4368469964806807, :y 1} {:x -0.5942926545724515, :y 3} {:x 1.4368469964806807, :y 1} {:x -0.4589393529282559, :y 2} {:x 1.4368469964806807, :y 1} {:x -0.4589393529282559, :y 2} {:x -0.4589393529282559, :y 2} {:x -0.5942926545724515, :y 3} {:x 1.4368469964806807, :y 1} {:x -0.4589393529282559, :y 2} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -0.5942926545724515, :y 3} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -0.5942926545724515, :y 3} {:x -0.8059581442333353, :y 4} {:x 1.4368469964806807, :y 1} {:x -0.8059581442333353, :y 4} {:x -0.4589393529282559, :y 2} {:x -0.4589393529282559, :y 2} {:x -0.5942926545724515, :y 3} {:x 1.4368469964806807, :y 1} {:x -0.4589393529282559, :y 2} {:x 1.4368469964806807, :y 1} {:x -0.8059581442333353, :y 4} {:x -0.5942926545724515, :y 3} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -0.4589393529282559, :y 2} {:x -0.5942926545724515, :y 3} {:x 1.4368469964806807, :y 1} {:x -0.4589393529282559, :y 2} {:x -0.8059581442333353, :y 4} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -0.4589393529282559, :y 2} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -0.4589393529282559, :y 2} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -0.5942926545724515, :y 3} {:x 1.4368469964806807, :y 1} {:x -0.4589393529282559, :y 2} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -0.8059581442333353, :y 4} {:x 1.4368469964806807, :y 1} {:x 1.4368469964806807, :y 1} {:x -2.0072115684121643, :y 8} {:x -2.0072115684121643, :y 8} {:x -2.0072115684121643, :y 8} {:x -2.0072115684121643, :y 8} {:x -1.9099417631646933, :y 7} {:x -1.9099417631646933, :y 7} {:x -1.9099417631646933, :y 7} {:x -1.9099417631646933, :y 7} {:x -1.7374705404343287, :y 6} {:x -1.7374705404343287, :y 6} {:x -1.7374705404343287, :y 6} {:x -1.7374705404343287, :y 6} {:x -1.452985831073355, :y 5} {:x -1.452985831073355, :y 5} {:x -1.452985831073355, :y 5} {:x -1.452985831073355, :y 5} {:x -2.298761224980861, :y 11} {:x -2.298761224980861, :y 11} {:x -2.17727991716631, :y 10} {:x -2.130653782248272, :y 9})}), :marks ({:type \"symbol\", :from {:data \"85f99bf7-f834-4b0f-bfab-4578e9502c38\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :fill {:value \"red\"}, :fillOpacity {:value 1}}, :update {:shape \"circle\", :size {:value 70}, :stroke {:value \"transparent\"}}, :hover {:size {:value 210}, :stroke {:value \"white\"}}}} {:type \"symbol\", :from {:data \"727d244f-8a9a-4fed-b01d-7404710b05ce\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :fill {:value \"blue\"}, :fillOpacity {:value 1}}, :update {:shape \"circle\", :size {:value 70}, :stroke {:value \"transparent\"}}, :hover {:size {:value 210}, :stroke {:value \"white\"}}}} {:type \"symbol\", :from {:data \"e3bcad74-936b-4b40-a31a-5fdf5d03a8c1\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :fill {:value \"#ff29d2\"}, :fillOpacity {:value 1}}, :update {:shape \"circle\", :size {:value 70}, :stroke {:value \"transparent\"}}, :hover {:size {:value 210}, :stroke {:value \"white\"}}}})}}"}
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
;;; {"type":"vega","content":{"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50},"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"17eaf8b0-5bf6-4852-b567-a9895dfb8c0e","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":[0,4705]}],"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"data":[{"name":"17eaf8b0-5bf6-4852-b567-a9895dfb8c0e","values":[{"x":0,"y":3031},{"x":1,"y":2637},{"x":2,"y":2838},{"x":3,"y":1869},{"x":4,"y":1351},{"x":5,"y":1734},{"x":6,"y":2059},{"x":7,"y":1612},{"x":8,"y":2100},{"x":9,"y":2165},{"x":10,"y":1197},{"x":11,"y":407},{"x":12,"y":959},{"x":13,"y":913},{"x":14,"y":848},{"x":15,"y":1009},{"x":16,"y":987},{"x":17,"y":890},{"x":18,"y":946},{"x":19,"y":966},{"x":20,"y":1400},{"x":21,"y":1718},{"x":22,"y":2031},{"x":23,"y":1650},{"x":24,"y":1348},{"x":25,"y":1258},{"x":26,"y":1490},{"x":27,"y":1899},{"x":28,"y":1674},{"x":29,"y":1100},{"x":30,"y":1399},{"x":31,"y":1214},{"x":32,"y":1460},{"x":33,"y":1084},{"x":34,"y":1215},{"x":35,"y":1594},{"x":36,"y":1194},{"x":37,"y":1239},{"x":38,"y":1771},{"x":39,"y":1526},{"x":40,"y":1051},{"x":41,"y":2470},{"x":42,"y":1664},{"x":43,"y":1225},{"x":44,"y":1762},{"x":45,"y":1188},{"x":46,"y":1796},{"x":47,"y":1497},{"x":48,"y":1914},{"x":49,"y":2222},{"x":50,"y":1302},{"x":51,"y":1136},{"x":52,"y":1460},{"x":53,"y":1611},{"x":54,"y":1331},{"x":55,"y":1318},{"x":56,"y":2424},{"x":57,"y":1350},{"x":58,"y":1283},{"x":59,"y":1921},{"x":60,"y":1702},{"x":61,"y":1019},{"x":62,"y":821},{"x":63,"y":1310},{"x":64,"y":1802},{"x":65,"y":2066},{"x":66,"y":1509},{"x":67,"y":2287},{"x":68,"y":2140},{"x":69,"y":1937},{"x":70,"y":1875},{"x":71,"y":1548},{"x":72,"y":1601},{"x":73,"y":2026},{"x":74,"y":1379},{"x":75,"y":1603},{"x":76,"y":2005},{"x":77,"y":1690},{"x":78,"y":1673},{"x":79,"y":1684},{"x":80,"y":1384},{"x":81,"y":1592},{"x":82,"y":1289},{"x":83,"y":920},{"x":84,"y":2628},{"x":85,"y":2592},{"x":86,"y":373},{"x":87,"y":545},{"x":88,"y":743},{"x":89,"y":1387},{"x":90,"y":744},{"x":91,"y":972},{"x":92,"y":1667},{"x":93,"y":2285},{"x":94,"y":1493},{"x":95,"y":1773},{"x":96,"y":1652},{"x":97,"y":821},{"x":98,"y":1127},{"x":99,"y":1863},{"x":100,"y":1325},{"x":101,"y":2153},{"x":102,"y":2383},{"x":103,"y":1900},{"x":104,"y":2430},{"x":105,"y":1021},{"x":106,"y":1964},{"x":107,"y":1894},{"x":108,"y":1609},{"x":109,"y":2735},{"x":110,"y":2225},{"x":111,"y":1434},{"x":112,"y":1407},{"x":113,"y":1866},{"x":114,"y":1911},{"x":115,"y":2302},{"x":116,"y":1728},{"x":117,"y":2156},{"x":118,"y":1484},{"x":119,"y":690},{"x":120,"y":1785},{"x":121,"y":1222},{"x":122,"y":1092},{"x":123,"y":881},{"x":124,"y":1377},{"x":125,"y":1120},{"x":126,"y":1207},{"x":127,"y":2382},{"x":128,"y":2127},{"x":129,"y":1442},{"x":130,"y":1421},{"x":131,"y":1703},{"x":132,"y":1918},{"x":133,"y":1344},{"x":134,"y":1465},{"x":135,"y":1353},{"x":136,"y":1304},{"x":137,"y":1654},{"x":138,"y":1873},{"x":139,"y":909},{"x":140,"y":1602},{"x":141,"y":1294},{"x":142,"y":1312},{"x":143,"y":1396},{"x":144,"y":1093},{"x":145,"y":1261},{"x":146,"y":1510},{"x":147,"y":1419},{"x":148,"y":1349},{"x":149,"y":1550},{"x":150,"y":1455},{"x":151,"y":1088},{"x":152,"y":759},{"x":153,"y":1376},{"x":154,"y":1554},{"x":155,"y":2036},{"x":156,"y":2141},{"x":157,"y":1327},{"x":158,"y":1599},{"x":159,"y":1744},{"x":160,"y":1822},{"x":161,"y":1973},{"x":162,"y":891},{"x":163,"y":1362},{"x":164,"y":1474},{"x":165,"y":1174},{"x":166,"y":1042},{"x":167,"y":992},{"x":168,"y":1336},{"x":169,"y":1381},{"x":170,"y":1492},{"x":171,"y":1796},{"x":172,"y":2131},{"x":173,"y":2566},{"x":174,"y":3189},{"x":175,"y":3030},{"x":176,"y":2253},{"x":177,"y":2453},{"x":178,"y":4705},{"x":179,"y":2415},{"x":180,"y":2477},{"x":181,"y":908},{"x":182,"y":1591},{"x":183,"y":3168},{"x":184,"y":1975},{"x":185,"y":1733},{"x":186,"y":1468},{"x":187,"y":2343},{"x":188,"y":4442},{"x":189,"y":2022},{"x":190,"y":1692},{"x":191,"y":2372},{"x":192,"y":1459},{"x":193,"y":1769},{"x":194,"y":1731},{"x":195,"y":1422},{"x":196,"y":2052},{"x":197,"y":1941},{"x":198,"y":2136},{"x":199,"y":1487}]},{"name":"d76170a5-22d3-412a-800a-608a1e3c2816","values":[{"x":0,"y":86},{"x":1,"y":147},{"x":2,"y":71},{"x":3,"y":78},{"x":4,"y":54},{"x":5,"y":53},{"x":6,"y":55},{"x":7,"y":49},{"x":8,"y":61},{"x":9,"y":49},{"x":10,"y":58},{"x":11,"y":56},{"x":12,"y":298},{"x":13,"y":365},{"x":14,"y":253},{"x":15,"y":305},{"x":16,"y":296},{"x":17,"y":290},{"x":18,"y":303},{"x":19,"y":243},{"x":20,"y":299},{"x":21,"y":207},{"x":22,"y":245},{"x":23,"y":228},{"x":24,"y":181},{"x":25,"y":197},{"x":26,"y":255},{"x":27,"y":193},{"x":28,"y":202},{"x":29,"y":222},{"x":30,"y":168},{"x":31,"y":188},{"x":32,"y":202},{"x":33,"y":162},{"x":34,"y":271},{"x":35,"y":260},{"x":36,"y":189},{"x":37,"y":157},{"x":38,"y":185},{"x":39,"y":212},{"x":40,"y":198},{"x":41,"y":240},{"x":42,"y":197},{"x":43,"y":164},{"x":44,"y":196},{"x":45,"y":215},{"x":46,"y":241},{"x":47,"y":192},{"x":48,"y":154},{"x":49,"y":173},{"x":50,"y":153},{"x":51,"y":240},{"x":52,"y":219},{"x":53,"y":297},{"x":54,"y":235},{"x":55,"y":180},{"x":56,"y":232},{"x":57,"y":120},{"x":58,"y":172},{"x":59,"y":233},{"x":60,"y":280},{"x":61,"y":246},{"x":62,"y":232},{"x":63,"y":193},{"x":64,"y":222},{"x":65,"y":227},{"x":66,"y":205},{"x":67,"y":171},{"x":68,"y":213},{"x":69,"y":206},{"x":70,"y":149},{"x":71,"y":205},{"x":72,"y":192},{"x":73,"y":234},{"x":74,"y":127},{"x":75,"y":228},{"x":76,"y":208},{"x":77,"y":125},{"x":78,"y":208},{"x":79,"y":183},{"x":80,"y":198},{"x":81,"y":125},{"x":82,"y":186},{"x":83,"y":229},{"x":84,"y":262},{"x":85,"y":206},{"x":86,"y":107},{"x":87,"y":267},{"x":88,"y":316},{"x":89,"y":258},{"x":90,"y":163},{"x":91,"y":218},{"x":92,"y":242},{"x":93,"y":275},{"x":94,"y":209},{"x":95,"y":217},{"x":96,"y":202},{"x":97,"y":237},{"x":98,"y":284},{"x":99,"y":193},{"x":100,"y":218},{"x":101,"y":168},{"x":102,"y":212},{"x":103,"y":128},{"x":104,"y":308},{"x":105,"y":237},{"x":106,"y":315},{"x":107,"y":258},{"x":108,"y":222},{"x":109,"y":229},{"x":110,"y":167},{"x":111,"y":159},{"x":112,"y":272},{"x":113,"y":309},{"x":114,"y":244},{"x":115,"y":207},{"x":116,"y":296},{"x":117,"y":327},{"x":118,"y":267},{"x":119,"y":259},{"x":120,"y":291},{"x":121,"y":286},{"x":122,"y":296},{"x":123,"y":380},{"x":124,"y":421},{"x":125,"y":320},{"x":126,"y":310},{"x":127,"y":271},{"x":128,"y":274},{"x":129,"y":201},{"x":130,"y":221},{"x":131,"y":208},{"x":132,"y":154},{"x":133,"y":168},{"x":134,"y":198},{"x":135,"y":154},{"x":136,"y":215},{"x":137,"y":273},{"x":138,"y":180},{"x":139,"y":119},{"x":140,"y":208},{"x":141,"y":212},{"x":142,"y":201},{"x":143,"y":206},{"x":144,"y":171},{"x":145,"y":204},{"x":146,"y":172},{"x":147,"y":206},{"x":148,"y":162},{"x":149,"y":140},{"x":150,"y":108},{"x":151,"y":167},{"x":152,"y":183},{"x":153,"y":234},{"x":154,"y":230},{"x":155,"y":179},{"x":156,"y":143},{"x":157,"y":151},{"x":158,"y":201},{"x":159,"y":166},{"x":160,"y":216},{"x":161,"y":236},{"x":162,"y":130},{"x":163,"y":215},{"x":164,"y":214},{"x":165,"y":164},{"x":166,"y":154},{"x":167,"y":206},{"x":168,"y":250},{"x":169,"y":227},{"x":170,"y":131},{"x":171,"y":163},{"x":172,"y":154},{"x":173,"y":250},{"x":174,"y":205},{"x":175,"y":217},{"x":176,"y":223},{"x":177,"y":227},{"x":178,"y":529},{"x":179,"y":206},{"x":180,"y":224},{"x":181,"y":243},{"x":182,"y":305},{"x":183,"y":214},{"x":184,"y":246},{"x":185,"y":198},{"x":186,"y":175},{"x":187,"y":179},{"x":188,"y":232},{"x":189,"y":140},{"x":190,"y":221},{"x":191,"y":190},{"x":192,"y":164},{"x":193,"y":148},{"x":194,"y":157},{"x":195,"y":203},{"x":196,"y":251},{"x":197,"y":149},{"x":198,"y":131},{"x":199,"y":183}]},{"name":"a2d08ed6-ffbe-496f-b540-c6142be30f24","values":[{"x":0,"y":25},{"x":1,"y":19},{"x":2,"y":4},{"x":3,"y":6},{"x":4,"y":6},{"x":5,"y":4},{"x":6,"y":3},{"x":7,"y":3},{"x":8,"y":3},{"x":9,"y":3},{"x":10,"y":4},{"x":11,"y":2},{"x":12,"y":3},{"x":13,"y":3},{"x":14,"y":2},{"x":15,"y":2},{"x":16,"y":2},{"x":17,"y":2},{"x":18,"y":1},{"x":19,"y":2},{"x":20,"y":3},{"x":21,"y":1},{"x":22,"y":2},{"x":23,"y":1},{"x":24,"y":1},{"x":25,"y":1},{"x":26,"y":1},{"x":27,"y":2},{"x":28,"y":2},{"x":29,"y":1},{"x":30,"y":1},{"x":31,"y":1},{"x":32,"y":2},{"x":33,"y":1},{"x":34,"y":1},{"x":35,"y":1},{"x":36,"y":0},{"x":37,"y":1},{"x":38,"y":1},{"x":39,"y":1},{"x":40,"y":1},{"x":41,"y":1},{"x":42,"y":1},{"x":43,"y":1},{"x":44,"y":1},{"x":45,"y":0},{"x":46,"y":1},{"x":47,"y":1},{"x":48,"y":2},{"x":49,"y":1},{"x":50,"y":0},{"x":51,"y":1},{"x":52,"y":1},{"x":53,"y":1},{"x":54,"y":1},{"x":55,"y":1},{"x":56,"y":1},{"x":57,"y":2},{"x":58,"y":1},{"x":59,"y":1},{"x":60,"y":1},{"x":61,"y":1},{"x":62,"y":1},{"x":63,"y":1},{"x":64,"y":2},{"x":65,"y":1},{"x":66,"y":1},{"x":67,"y":1},{"x":68,"y":1},{"x":69,"y":1},{"x":70,"y":2},{"x":71,"y":1},{"x":72,"y":1},{"x":73,"y":1},{"x":74,"y":0},{"x":75,"y":1},{"x":76,"y":0},{"x":77,"y":1},{"x":78,"y":1},{"x":79,"y":1},{"x":80,"y":1},{"x":81,"y":1},{"x":82,"y":1},{"x":83,"y":0},{"x":84,"y":1},{"x":85,"y":1},{"x":86,"y":0},{"x":87,"y":1},{"x":88,"y":0},{"x":89,"y":1},{"x":90,"y":1},{"x":91,"y":1},{"x":92,"y":2},{"x":93,"y":1},{"x":94,"y":0},{"x":95,"y":1},{"x":96,"y":1},{"x":97,"y":1},{"x":98,"y":1},{"x":99,"y":1},{"x":100,"y":0},{"x":101,"y":1},{"x":102,"y":1},{"x":103,"y":1},{"x":104,"y":1},{"x":105,"y":1},{"x":106,"y":0},{"x":107,"y":2},{"x":108,"y":1},{"x":109,"y":1},{"x":110,"y":1},{"x":111,"y":2},{"x":112,"y":1},{"x":113,"y":1},{"x":114,"y":1},{"x":115,"y":1},{"x":116,"y":1},{"x":117,"y":1},{"x":118,"y":1},{"x":119,"y":1},{"x":120,"y":1},{"x":121,"y":1},{"x":122,"y":1},{"x":123,"y":1},{"x":124,"y":1},{"x":125,"y":1},{"x":126,"y":0},{"x":127,"y":1},{"x":128,"y":1},{"x":129,"y":1},{"x":130,"y":0},{"x":131,"y":1},{"x":132,"y":1},{"x":133,"y":1},{"x":134,"y":1},{"x":135,"y":1},{"x":136,"y":0},{"x":137,"y":1},{"x":138,"y":0},{"x":139,"y":0},{"x":140,"y":0},{"x":141,"y":1},{"x":142,"y":1},{"x":143,"y":2},{"x":144,"y":1},{"x":145,"y":0},{"x":146,"y":1},{"x":147,"y":1},{"x":148,"y":1},{"x":149,"y":1},{"x":150,"y":1},{"x":151,"y":0},{"x":152,"y":1},{"x":153,"y":1},{"x":154,"y":1},{"x":155,"y":0},{"x":156,"y":1},{"x":157,"y":1},{"x":158,"y":1},{"x":159,"y":1},{"x":160,"y":1},{"x":161,"y":1},{"x":162,"y":1},{"x":163,"y":1},{"x":164,"y":0},{"x":165,"y":1},{"x":166,"y":0},{"x":167,"y":0},{"x":168,"y":0},{"x":169,"y":0},{"x":170,"y":1},{"x":171,"y":1},{"x":172,"y":1},{"x":173,"y":1},{"x":174,"y":1},{"x":175,"y":1},{"x":176,"y":2},{"x":177,"y":1},{"x":178,"y":12},{"x":179,"y":1},{"x":180,"y":5},{"x":181,"y":1},{"x":182,"y":0},{"x":183,"y":2},{"x":184,"y":4},{"x":185,"y":1},{"x":186,"y":1},{"x":187,"y":2},{"x":188,"y":1},{"x":189,"y":1},{"x":190,"y":2},{"x":191,"y":1},{"x":192,"y":1},{"x":193,"y":1},{"x":194,"y":0},{"x":195,"y":1},{"x":196,"y":1},{"x":197,"y":1},{"x":198,"y":1},{"x":199,"y":1}]},{"name":"8c3d5841-36d8-4dc0-bf2e-1bf47a7de4e1","values":[{"x":0,"y":2920},{"x":1,"y":2471},{"x":2,"y":2763},{"x":3,"y":1785},{"x":4,"y":1291},{"x":5,"y":1677},{"x":6,"y":2001},{"x":7,"y":1560},{"x":8,"y":2036},{"x":9,"y":2113},{"x":10,"y":1135},{"x":11,"y":349},{"x":12,"y":658},{"x":13,"y":545},{"x":14,"y":593},{"x":15,"y":702},{"x":16,"y":689},{"x":17,"y":598},{"x":18,"y":642},{"x":19,"y":721},{"x":20,"y":1098},{"x":21,"y":1510},{"x":22,"y":1784},{"x":23,"y":1421},{"x":24,"y":1166},{"x":25,"y":1060},{"x":26,"y":1234},{"x":27,"y":1704},{"x":28,"y":1470},{"x":29,"y":877},{"x":30,"y":1230},{"x":31,"y":1025},{"x":32,"y":1256},{"x":33,"y":921},{"x":34,"y":943},{"x":35,"y":1333},{"x":36,"y":1005},{"x":37,"y":1081},{"x":38,"y":1585},{"x":39,"y":1313},{"x":40,"y":852},{"x":41,"y":2229},{"x":42,"y":1466},{"x":43,"y":1060},{"x":44,"y":1565},{"x":45,"y":973},{"x":46,"y":1554},{"x":47,"y":1304},{"x":48,"y":1758},{"x":49,"y":2048},{"x":50,"y":1149},{"x":51,"y":895},{"x":52,"y":1240},{"x":53,"y":1313},{"x":54,"y":1095},{"x":55,"y":1137},{"x":56,"y":2191},{"x":57,"y":1228},{"x":58,"y":1110},{"x":59,"y":1687},{"x":60,"y":1421},{"x":61,"y":772},{"x":62,"y":588},{"x":63,"y":1116},{"x":64,"y":1578},{"x":65,"y":1838},{"x":66,"y":1303},{"x":67,"y":2115},{"x":68,"y":1926},{"x":69,"y":1730},{"x":70,"y":1724},{"x":71,"y":1342},{"x":72,"y":1408},{"x":73,"y":1791},{"x":74,"y":1252},{"x":75,"y":1374},{"x":76,"y":1797},{"x":77,"y":1564},{"x":78,"y":1464},{"x":79,"y":1500},{"x":80,"y":1185},{"x":81,"y":1466},{"x":82,"y":1102},{"x":83,"y":691},{"x":84,"y":2365},{"x":85,"y":2385},{"x":86,"y":266},{"x":87,"y":277},{"x":88,"y":427},{"x":89,"y":1128},{"x":90,"y":580},{"x":91,"y":753},{"x":92,"y":1423},{"x":93,"y":2009},{"x":94,"y":1284},{"x":95,"y":1555},{"x":96,"y":1449},{"x":97,"y":583},{"x":98,"y":842},{"x":99,"y":1669},{"x":100,"y":1107},{"x":101,"y":1984},{"x":102,"y":2170},{"x":103,"y":1771},{"x":104,"y":2121},{"x":105,"y":783},{"x":106,"y":1649},{"x":107,"y":1634},{"x":108,"y":1386},{"x":109,"y":2505},{"x":110,"y":2057},{"x":111,"y":1273},{"x":112,"y":1134},{"x":113,"y":1556},{"x":114,"y":1666},{"x":115,"y":2094},{"x":116,"y":1431},{"x":117,"y":1828},{"x":118,"y":1216},{"x":119,"y":430},{"x":120,"y":1493},{"x":121,"y":935},{"x":122,"y":795},{"x":123,"y":500},{"x":124,"y":955},{"x":125,"y":799},{"x":126,"y":897},{"x":127,"y":2110},{"x":128,"y":1852},{"x":129,"y":1240},{"x":130,"y":1200},{"x":131,"y":1494},{"x":132,"y":1763},{"x":133,"y":1175},{"x":134,"y":1266},{"x":135,"y":1198},{"x":136,"y":1089},{"x":137,"y":1380},{"x":138,"y":1693},{"x":139,"y":790},{"x":140,"y":1394},{"x":141,"y":1081},{"x":142,"y":1110},{"x":143,"y":1188},{"x":144,"y":921},{"x":145,"y":1057},{"x":146,"y":1337},{"x":147,"y":1212},{"x":148,"y":1186},{"x":149,"y":1409},{"x":150,"y":1346},{"x":151,"y":921},{"x":152,"y":575},{"x":153,"y":1141},{"x":154,"y":1323},{"x":155,"y":1857},{"x":156,"y":1997},{"x":157,"y":1175},{"x":158,"y":1397},{"x":159,"y":1577},{"x":160,"y":1605},{"x":161,"y":1736},{"x":162,"y":760},{"x":163,"y":1146},{"x":164,"y":1260},{"x":165,"y":1009},{"x":166,"y":888},{"x":167,"y":786},{"x":168,"y":1086},{"x":169,"y":1154},{"x":170,"y":1360},{"x":171,"y":1632},{"x":172,"y":1976},{"x":173,"y":2315},{"x":174,"y":2983},{"x":175,"y":2812},{"x":176,"y":2028},{"x":177,"y":2225},{"x":178,"y":4164},{"x":179,"y":2208},{"x":180,"y":2248},{"x":181,"y":664},{"x":182,"y":1286},{"x":183,"y":2952},{"x":184,"y":1725},{"x":185,"y":1534},{"x":186,"y":1292},{"x":187,"y":2162},{"x":188,"y":4209},{"x":189,"y":1881},{"x":190,"y":1469},{"x":191,"y":2181},{"x":192,"y":1294},{"x":193,"y":1620},{"x":194,"y":1574},{"x":195,"y":1218},{"x":196,"y":1800},{"x":197,"y":1791},{"x":198,"y":2004},{"x":199,"y":1303}]}],"marks":[{"type":"line","from":{"data":"17eaf8b0-5bf6-4852-b567-a9895dfb8c0e"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"yellow"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}},{"type":"line","from":{"data":"d76170a5-22d3-412a-800a-608a1e3c2816"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"red"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}},{"type":"line","from":{"data":"a2d08ed6-ffbe-496f-b540-c6142be30f24"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"green"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}},{"type":"line","from":{"data":"8c3d5841-36d8-4dc0-bf2e-1bf47a7de4e1"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"blue"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}}]},"value":"#gorilla_repl.vega.VegaView{:content {:width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}, :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"17eaf8b0-5bf6-4852-b567-a9895dfb8c0e\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain [0 4705]}], :axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :data ({:name \"17eaf8b0-5bf6-4852-b567-a9895dfb8c0e\", :values ({:x 0, :y 3031} {:x 1, :y 2637} {:x 2, :y 2838} {:x 3, :y 1869} {:x 4, :y 1351} {:x 5, :y 1734} {:x 6, :y 2059} {:x 7, :y 1612} {:x 8, :y 2100} {:x 9, :y 2165} {:x 10, :y 1197} {:x 11, :y 407} {:x 12, :y 959} {:x 13, :y 913} {:x 14, :y 848} {:x 15, :y 1009} {:x 16, :y 987} {:x 17, :y 890} {:x 18, :y 946} {:x 19, :y 966} {:x 20, :y 1400} {:x 21, :y 1718} {:x 22, :y 2031} {:x 23, :y 1650} {:x 24, :y 1348} {:x 25, :y 1258} {:x 26, :y 1490} {:x 27, :y 1899} {:x 28, :y 1674} {:x 29, :y 1100} {:x 30, :y 1399} {:x 31, :y 1214} {:x 32, :y 1460} {:x 33, :y 1084} {:x 34, :y 1215} {:x 35, :y 1594} {:x 36, :y 1194} {:x 37, :y 1239} {:x 38, :y 1771} {:x 39, :y 1526} {:x 40, :y 1051} {:x 41, :y 2470} {:x 42, :y 1664} {:x 43, :y 1225} {:x 44, :y 1762} {:x 45, :y 1188} {:x 46, :y 1796} {:x 47, :y 1497} {:x 48, :y 1914} {:x 49, :y 2222} {:x 50, :y 1302} {:x 51, :y 1136} {:x 52, :y 1460} {:x 53, :y 1611} {:x 54, :y 1331} {:x 55, :y 1318} {:x 56, :y 2424} {:x 57, :y 1350} {:x 58, :y 1283} {:x 59, :y 1921} {:x 60, :y 1702} {:x 61, :y 1019} {:x 62, :y 821} {:x 63, :y 1310} {:x 64, :y 1802} {:x 65, :y 2066} {:x 66, :y 1509} {:x 67, :y 2287} {:x 68, :y 2140} {:x 69, :y 1937} {:x 70, :y 1875} {:x 71, :y 1548} {:x 72, :y 1601} {:x 73, :y 2026} {:x 74, :y 1379} {:x 75, :y 1603} {:x 76, :y 2005} {:x 77, :y 1690} {:x 78, :y 1673} {:x 79, :y 1684} {:x 80, :y 1384} {:x 81, :y 1592} {:x 82, :y 1289} {:x 83, :y 920} {:x 84, :y 2628} {:x 85, :y 2592} {:x 86, :y 373} {:x 87, :y 545} {:x 88, :y 743} {:x 89, :y 1387} {:x 90, :y 744} {:x 91, :y 972} {:x 92, :y 1667} {:x 93, :y 2285} {:x 94, :y 1493} {:x 95, :y 1773} {:x 96, :y 1652} {:x 97, :y 821} {:x 98, :y 1127} {:x 99, :y 1863} {:x 100, :y 1325} {:x 101, :y 2153} {:x 102, :y 2383} {:x 103, :y 1900} {:x 104, :y 2430} {:x 105, :y 1021} {:x 106, :y 1964} {:x 107, :y 1894} {:x 108, :y 1609} {:x 109, :y 2735} {:x 110, :y 2225} {:x 111, :y 1434} {:x 112, :y 1407} {:x 113, :y 1866} {:x 114, :y 1911} {:x 115, :y 2302} {:x 116, :y 1728} {:x 117, :y 2156} {:x 118, :y 1484} {:x 119, :y 690} {:x 120, :y 1785} {:x 121, :y 1222} {:x 122, :y 1092} {:x 123, :y 881} {:x 124, :y 1377} {:x 125, :y 1120} {:x 126, :y 1207} {:x 127, :y 2382} {:x 128, :y 2127} {:x 129, :y 1442} {:x 130, :y 1421} {:x 131, :y 1703} {:x 132, :y 1918} {:x 133, :y 1344} {:x 134, :y 1465} {:x 135, :y 1353} {:x 136, :y 1304} {:x 137, :y 1654} {:x 138, :y 1873} {:x 139, :y 909} {:x 140, :y 1602} {:x 141, :y 1294} {:x 142, :y 1312} {:x 143, :y 1396} {:x 144, :y 1093} {:x 145, :y 1261} {:x 146, :y 1510} {:x 147, :y 1419} {:x 148, :y 1349} {:x 149, :y 1550} {:x 150, :y 1455} {:x 151, :y 1088} {:x 152, :y 759} {:x 153, :y 1376} {:x 154, :y 1554} {:x 155, :y 2036} {:x 156, :y 2141} {:x 157, :y 1327} {:x 158, :y 1599} {:x 159, :y 1744} {:x 160, :y 1822} {:x 161, :y 1973} {:x 162, :y 891} {:x 163, :y 1362} {:x 164, :y 1474} {:x 165, :y 1174} {:x 166, :y 1042} {:x 167, :y 992} {:x 168, :y 1336} {:x 169, :y 1381} {:x 170, :y 1492} {:x 171, :y 1796} {:x 172, :y 2131} {:x 173, :y 2566} {:x 174, :y 3189} {:x 175, :y 3030} {:x 176, :y 2253} {:x 177, :y 2453} {:x 178, :y 4705} {:x 179, :y 2415} {:x 180, :y 2477} {:x 181, :y 908} {:x 182, :y 1591} {:x 183, :y 3168} {:x 184, :y 1975} {:x 185, :y 1733} {:x 186, :y 1468} {:x 187, :y 2343} {:x 188, :y 4442} {:x 189, :y 2022} {:x 190, :y 1692} {:x 191, :y 2372} {:x 192, :y 1459} {:x 193, :y 1769} {:x 194, :y 1731} {:x 195, :y 1422} {:x 196, :y 2052} {:x 197, :y 1941} {:x 198, :y 2136} {:x 199, :y 1487})} {:name \"d76170a5-22d3-412a-800a-608a1e3c2816\", :values ({:x 0, :y 86} {:x 1, :y 147} {:x 2, :y 71} {:x 3, :y 78} {:x 4, :y 54} {:x 5, :y 53} {:x 6, :y 55} {:x 7, :y 49} {:x 8, :y 61} {:x 9, :y 49} {:x 10, :y 58} {:x 11, :y 56} {:x 12, :y 298} {:x 13, :y 365} {:x 14, :y 253} {:x 15, :y 305} {:x 16, :y 296} {:x 17, :y 290} {:x 18, :y 303} {:x 19, :y 243} {:x 20, :y 299} {:x 21, :y 207} {:x 22, :y 245} {:x 23, :y 228} {:x 24, :y 181} {:x 25, :y 197} {:x 26, :y 255} {:x 27, :y 193} {:x 28, :y 202} {:x 29, :y 222} {:x 30, :y 168} {:x 31, :y 188} {:x 32, :y 202} {:x 33, :y 162} {:x 34, :y 271} {:x 35, :y 260} {:x 36, :y 189} {:x 37, :y 157} {:x 38, :y 185} {:x 39, :y 212} {:x 40, :y 198} {:x 41, :y 240} {:x 42, :y 197} {:x 43, :y 164} {:x 44, :y 196} {:x 45, :y 215} {:x 46, :y 241} {:x 47, :y 192} {:x 48, :y 154} {:x 49, :y 173} {:x 50, :y 153} {:x 51, :y 240} {:x 52, :y 219} {:x 53, :y 297} {:x 54, :y 235} {:x 55, :y 180} {:x 56, :y 232} {:x 57, :y 120} {:x 58, :y 172} {:x 59, :y 233} {:x 60, :y 280} {:x 61, :y 246} {:x 62, :y 232} {:x 63, :y 193} {:x 64, :y 222} {:x 65, :y 227} {:x 66, :y 205} {:x 67, :y 171} {:x 68, :y 213} {:x 69, :y 206} {:x 70, :y 149} {:x 71, :y 205} {:x 72, :y 192} {:x 73, :y 234} {:x 74, :y 127} {:x 75, :y 228} {:x 76, :y 208} {:x 77, :y 125} {:x 78, :y 208} {:x 79, :y 183} {:x 80, :y 198} {:x 81, :y 125} {:x 82, :y 186} {:x 83, :y 229} {:x 84, :y 262} {:x 85, :y 206} {:x 86, :y 107} {:x 87, :y 267} {:x 88, :y 316} {:x 89, :y 258} {:x 90, :y 163} {:x 91, :y 218} {:x 92, :y 242} {:x 93, :y 275} {:x 94, :y 209} {:x 95, :y 217} {:x 96, :y 202} {:x 97, :y 237} {:x 98, :y 284} {:x 99, :y 193} {:x 100, :y 218} {:x 101, :y 168} {:x 102, :y 212} {:x 103, :y 128} {:x 104, :y 308} {:x 105, :y 237} {:x 106, :y 315} {:x 107, :y 258} {:x 108, :y 222} {:x 109, :y 229} {:x 110, :y 167} {:x 111, :y 159} {:x 112, :y 272} {:x 113, :y 309} {:x 114, :y 244} {:x 115, :y 207} {:x 116, :y 296} {:x 117, :y 327} {:x 118, :y 267} {:x 119, :y 259} {:x 120, :y 291} {:x 121, :y 286} {:x 122, :y 296} {:x 123, :y 380} {:x 124, :y 421} {:x 125, :y 320} {:x 126, :y 310} {:x 127, :y 271} {:x 128, :y 274} {:x 129, :y 201} {:x 130, :y 221} {:x 131, :y 208} {:x 132, :y 154} {:x 133, :y 168} {:x 134, :y 198} {:x 135, :y 154} {:x 136, :y 215} {:x 137, :y 273} {:x 138, :y 180} {:x 139, :y 119} {:x 140, :y 208} {:x 141, :y 212} {:x 142, :y 201} {:x 143, :y 206} {:x 144, :y 171} {:x 145, :y 204} {:x 146, :y 172} {:x 147, :y 206} {:x 148, :y 162} {:x 149, :y 140} {:x 150, :y 108} {:x 151, :y 167} {:x 152, :y 183} {:x 153, :y 234} {:x 154, :y 230} {:x 155, :y 179} {:x 156, :y 143} {:x 157, :y 151} {:x 158, :y 201} {:x 159, :y 166} {:x 160, :y 216} {:x 161, :y 236} {:x 162, :y 130} {:x 163, :y 215} {:x 164, :y 214} {:x 165, :y 164} {:x 166, :y 154} {:x 167, :y 206} {:x 168, :y 250} {:x 169, :y 227} {:x 170, :y 131} {:x 171, :y 163} {:x 172, :y 154} {:x 173, :y 250} {:x 174, :y 205} {:x 175, :y 217} {:x 176, :y 223} {:x 177, :y 227} {:x 178, :y 529} {:x 179, :y 206} {:x 180, :y 224} {:x 181, :y 243} {:x 182, :y 305} {:x 183, :y 214} {:x 184, :y 246} {:x 185, :y 198} {:x 186, :y 175} {:x 187, :y 179} {:x 188, :y 232} {:x 189, :y 140} {:x 190, :y 221} {:x 191, :y 190} {:x 192, :y 164} {:x 193, :y 148} {:x 194, :y 157} {:x 195, :y 203} {:x 196, :y 251} {:x 197, :y 149} {:x 198, :y 131} {:x 199, :y 183})} {:name \"a2d08ed6-ffbe-496f-b540-c6142be30f24\", :values ({:x 0, :y 25} {:x 1, :y 19} {:x 2, :y 4} {:x 3, :y 6} {:x 4, :y 6} {:x 5, :y 4} {:x 6, :y 3} {:x 7, :y 3} {:x 8, :y 3} {:x 9, :y 3} {:x 10, :y 4} {:x 11, :y 2} {:x 12, :y 3} {:x 13, :y 3} {:x 14, :y 2} {:x 15, :y 2} {:x 16, :y 2} {:x 17, :y 2} {:x 18, :y 1} {:x 19, :y 2} {:x 20, :y 3} {:x 21, :y 1} {:x 22, :y 2} {:x 23, :y 1} {:x 24, :y 1} {:x 25, :y 1} {:x 26, :y 1} {:x 27, :y 2} {:x 28, :y 2} {:x 29, :y 1} {:x 30, :y 1} {:x 31, :y 1} {:x 32, :y 2} {:x 33, :y 1} {:x 34, :y 1} {:x 35, :y 1} {:x 36, :y 0} {:x 37, :y 1} {:x 38, :y 1} {:x 39, :y 1} {:x 40, :y 1} {:x 41, :y 1} {:x 42, :y 1} {:x 43, :y 1} {:x 44, :y 1} {:x 45, :y 0} {:x 46, :y 1} {:x 47, :y 1} {:x 48, :y 2} {:x 49, :y 1} {:x 50, :y 0} {:x 51, :y 1} {:x 52, :y 1} {:x 53, :y 1} {:x 54, :y 1} {:x 55, :y 1} {:x 56, :y 1} {:x 57, :y 2} {:x 58, :y 1} {:x 59, :y 1} {:x 60, :y 1} {:x 61, :y 1} {:x 62, :y 1} {:x 63, :y 1} {:x 64, :y 2} {:x 65, :y 1} {:x 66, :y 1} {:x 67, :y 1} {:x 68, :y 1} {:x 69, :y 1} {:x 70, :y 2} {:x 71, :y 1} {:x 72, :y 1} {:x 73, :y 1} {:x 74, :y 0} {:x 75, :y 1} {:x 76, :y 0} {:x 77, :y 1} {:x 78, :y 1} {:x 79, :y 1} {:x 80, :y 1} {:x 81, :y 1} {:x 82, :y 1} {:x 83, :y 0} {:x 84, :y 1} {:x 85, :y 1} {:x 86, :y 0} {:x 87, :y 1} {:x 88, :y 0} {:x 89, :y 1} {:x 90, :y 1} {:x 91, :y 1} {:x 92, :y 2} {:x 93, :y 1} {:x 94, :y 0} {:x 95, :y 1} {:x 96, :y 1} {:x 97, :y 1} {:x 98, :y 1} {:x 99, :y 1} {:x 100, :y 0} {:x 101, :y 1} {:x 102, :y 1} {:x 103, :y 1} {:x 104, :y 1} {:x 105, :y 1} {:x 106, :y 0} {:x 107, :y 2} {:x 108, :y 1} {:x 109, :y 1} {:x 110, :y 1} {:x 111, :y 2} {:x 112, :y 1} {:x 113, :y 1} {:x 114, :y 1} {:x 115, :y 1} {:x 116, :y 1} {:x 117, :y 1} {:x 118, :y 1} {:x 119, :y 1} {:x 120, :y 1} {:x 121, :y 1} {:x 122, :y 1} {:x 123, :y 1} {:x 124, :y 1} {:x 125, :y 1} {:x 126, :y 0} {:x 127, :y 1} {:x 128, :y 1} {:x 129, :y 1} {:x 130, :y 0} {:x 131, :y 1} {:x 132, :y 1} {:x 133, :y 1} {:x 134, :y 1} {:x 135, :y 1} {:x 136, :y 0} {:x 137, :y 1} {:x 138, :y 0} {:x 139, :y 0} {:x 140, :y 0} {:x 141, :y 1} {:x 142, :y 1} {:x 143, :y 2} {:x 144, :y 1} {:x 145, :y 0} {:x 146, :y 1} {:x 147, :y 1} {:x 148, :y 1} {:x 149, :y 1} {:x 150, :y 1} {:x 151, :y 0} {:x 152, :y 1} {:x 153, :y 1} {:x 154, :y 1} {:x 155, :y 0} {:x 156, :y 1} {:x 157, :y 1} {:x 158, :y 1} {:x 159, :y 1} {:x 160, :y 1} {:x 161, :y 1} {:x 162, :y 1} {:x 163, :y 1} {:x 164, :y 0} {:x 165, :y 1} {:x 166, :y 0} {:x 167, :y 0} {:x 168, :y 0} {:x 169, :y 0} {:x 170, :y 1} {:x 171, :y 1} {:x 172, :y 1} {:x 173, :y 1} {:x 174, :y 1} {:x 175, :y 1} {:x 176, :y 2} {:x 177, :y 1} {:x 178, :y 12} {:x 179, :y 1} {:x 180, :y 5} {:x 181, :y 1} {:x 182, :y 0} {:x 183, :y 2} {:x 184, :y 4} {:x 185, :y 1} {:x 186, :y 1} {:x 187, :y 2} {:x 188, :y 1} {:x 189, :y 1} {:x 190, :y 2} {:x 191, :y 1} {:x 192, :y 1} {:x 193, :y 1} {:x 194, :y 0} {:x 195, :y 1} {:x 196, :y 1} {:x 197, :y 1} {:x 198, :y 1} {:x 199, :y 1})} {:name \"8c3d5841-36d8-4dc0-bf2e-1bf47a7de4e1\", :values ({:x 0, :y 2920} {:x 1, :y 2471} {:x 2, :y 2763} {:x 3, :y 1785} {:x 4, :y 1291} {:x 5, :y 1677} {:x 6, :y 2001} {:x 7, :y 1560} {:x 8, :y 2036} {:x 9, :y 2113} {:x 10, :y 1135} {:x 11, :y 349} {:x 12, :y 658} {:x 13, :y 545} {:x 14, :y 593} {:x 15, :y 702} {:x 16, :y 689} {:x 17, :y 598} {:x 18, :y 642} {:x 19, :y 721} {:x 20, :y 1098} {:x 21, :y 1510} {:x 22, :y 1784} {:x 23, :y 1421} {:x 24, :y 1166} {:x 25, :y 1060} {:x 26, :y 1234} {:x 27, :y 1704} {:x 28, :y 1470} {:x 29, :y 877} {:x 30, :y 1230} {:x 31, :y 1025} {:x 32, :y 1256} {:x 33, :y 921} {:x 34, :y 943} {:x 35, :y 1333} {:x 36, :y 1005} {:x 37, :y 1081} {:x 38, :y 1585} {:x 39, :y 1313} {:x 40, :y 852} {:x 41, :y 2229} {:x 42, :y 1466} {:x 43, :y 1060} {:x 44, :y 1565} {:x 45, :y 973} {:x 46, :y 1554} {:x 47, :y 1304} {:x 48, :y 1758} {:x 49, :y 2048} {:x 50, :y 1149} {:x 51, :y 895} {:x 52, :y 1240} {:x 53, :y 1313} {:x 54, :y 1095} {:x 55, :y 1137} {:x 56, :y 2191} {:x 57, :y 1228} {:x 58, :y 1110} {:x 59, :y 1687} {:x 60, :y 1421} {:x 61, :y 772} {:x 62, :y 588} {:x 63, :y 1116} {:x 64, :y 1578} {:x 65, :y 1838} {:x 66, :y 1303} {:x 67, :y 2115} {:x 68, :y 1926} {:x 69, :y 1730} {:x 70, :y 1724} {:x 71, :y 1342} {:x 72, :y 1408} {:x 73, :y 1791} {:x 74, :y 1252} {:x 75, :y 1374} {:x 76, :y 1797} {:x 77, :y 1564} {:x 78, :y 1464} {:x 79, :y 1500} {:x 80, :y 1185} {:x 81, :y 1466} {:x 82, :y 1102} {:x 83, :y 691} {:x 84, :y 2365} {:x 85, :y 2385} {:x 86, :y 266} {:x 87, :y 277} {:x 88, :y 427} {:x 89, :y 1128} {:x 90, :y 580} {:x 91, :y 753} {:x 92, :y 1423} {:x 93, :y 2009} {:x 94, :y 1284} {:x 95, :y 1555} {:x 96, :y 1449} {:x 97, :y 583} {:x 98, :y 842} {:x 99, :y 1669} {:x 100, :y 1107} {:x 101, :y 1984} {:x 102, :y 2170} {:x 103, :y 1771} {:x 104, :y 2121} {:x 105, :y 783} {:x 106, :y 1649} {:x 107, :y 1634} {:x 108, :y 1386} {:x 109, :y 2505} {:x 110, :y 2057} {:x 111, :y 1273} {:x 112, :y 1134} {:x 113, :y 1556} {:x 114, :y 1666} {:x 115, :y 2094} {:x 116, :y 1431} {:x 117, :y 1828} {:x 118, :y 1216} {:x 119, :y 430} {:x 120, :y 1493} {:x 121, :y 935} {:x 122, :y 795} {:x 123, :y 500} {:x 124, :y 955} {:x 125, :y 799} {:x 126, :y 897} {:x 127, :y 2110} {:x 128, :y 1852} {:x 129, :y 1240} {:x 130, :y 1200} {:x 131, :y 1494} {:x 132, :y 1763} {:x 133, :y 1175} {:x 134, :y 1266} {:x 135, :y 1198} {:x 136, :y 1089} {:x 137, :y 1380} {:x 138, :y 1693} {:x 139, :y 790} {:x 140, :y 1394} {:x 141, :y 1081} {:x 142, :y 1110} {:x 143, :y 1188} {:x 144, :y 921} {:x 145, :y 1057} {:x 146, :y 1337} {:x 147, :y 1212} {:x 148, :y 1186} {:x 149, :y 1409} {:x 150, :y 1346} {:x 151, :y 921} {:x 152, :y 575} {:x 153, :y 1141} {:x 154, :y 1323} {:x 155, :y 1857} {:x 156, :y 1997} {:x 157, :y 1175} {:x 158, :y 1397} {:x 159, :y 1577} {:x 160, :y 1605} {:x 161, :y 1736} {:x 162, :y 760} {:x 163, :y 1146} {:x 164, :y 1260} {:x 165, :y 1009} {:x 166, :y 888} {:x 167, :y 786} {:x 168, :y 1086} {:x 169, :y 1154} {:x 170, :y 1360} {:x 171, :y 1632} {:x 172, :y 1976} {:x 173, :y 2315} {:x 174, :y 2983} {:x 175, :y 2812} {:x 176, :y 2028} {:x 177, :y 2225} {:x 178, :y 4164} {:x 179, :y 2208} {:x 180, :y 2248} {:x 181, :y 664} {:x 182, :y 1286} {:x 183, :y 2952} {:x 184, :y 1725} {:x 185, :y 1534} {:x 186, :y 1292} {:x 187, :y 2162} {:x 188, :y 4209} {:x 189, :y 1881} {:x 190, :y 1469} {:x 191, :y 2181} {:x 192, :y 1294} {:x 193, :y 1620} {:x 194, :y 1574} {:x 195, :y 1218} {:x 196, :y 1800} {:x 197, :y 1791} {:x 198, :y 2004} {:x 199, :y 1303})}), :marks ({:type \"line\", :from {:data \"17eaf8b0-5bf6-4852-b567-a9895dfb8c0e\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"yellow\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}} {:type \"line\", :from {:data \"d76170a5-22d3-412a-800a-608a1e3c2816\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"red\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}} {:type \"line\", :from {:data \"a2d08ed6-ffbe-496f-b540-c6142be30f24\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"green\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}} {:type \"line\", :from {:data \"8c3d5841-36d8-4dc0-bf2e-1bf47a7de4e1\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"blue\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}})}}"}
;; <=

;; **
;;; Spread of complexity in population
;; **

;; @@
(plot/histogram (map :complexity (:rabble result)))
;; @@
;; =>
;;; {"type":"vega","content":{"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"f50f1304-d1f3-41cb-8a7f-1600db887930","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":{"data":"f50f1304-d1f3-41cb-8a7f-1600db887930","field":"data.y"}}],"marks":[{"type":"line","from":{"data":"f50f1304-d1f3-41cb-8a7f-1600db887930"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"interpolate":{"value":"step-before"},"fill":{"value":"steelblue"},"fillOpacity":{"value":0.4},"stroke":{"value":"steelblue"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}}],"data":[{"name":"f50f1304-d1f3-41cb-8a7f-1600db887930","values":[{"x":1.0,"y":0},{"x":2.125,"y":44.0},{"x":3.25,"y":9.0},{"x":4.375,"y":14.0},{"x":5.5,"y":8.0},{"x":6.625,"y":9.0},{"x":7.75,"y":3.0},{"x":8.875,"y":8.0},{"x":10.0,"y":2.0},{"x":11.125,"y":3.0},{"x":12.25,"y":0}]}],"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50}},"value":"#gorilla_repl.vega.VegaView{:content {:axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"f50f1304-d1f3-41cb-8a7f-1600db887930\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain {:data \"f50f1304-d1f3-41cb-8a7f-1600db887930\", :field \"data.y\"}}], :marks [{:type \"line\", :from {:data \"f50f1304-d1f3-41cb-8a7f-1600db887930\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :interpolate {:value \"step-before\"}, :fill {:value \"steelblue\"}, :fillOpacity {:value 0.4}, :stroke {:value \"steelblue\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}}], :data [{:name \"f50f1304-d1f3-41cb-8a7f-1600db887930\", :values ({:x 1.0, :y 0} {:x 2.125, :y 44.0} {:x 3.25, :y 9.0} {:x 4.375, :y 14.0} {:x 5.5, :y 8.0} {:x 6.625, :y 9.0} {:x 7.75, :y 3.0} {:x 8.875, :y 8.0} {:x 10.0, :y 2.0} {:x 11.125, :y 3.0} {:x 12.25, :y 0})}], :width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}}}"}
;; <=

;; @@
(plot/list-plot (:mean (:complexity @metrics/metrics)) :joined true)
;; @@
;; =>
;;; {"type":"vega","content":{"axes":[{"scale":"x","type":"x"},{"scale":"y","type":"y"}],"scales":[{"name":"x","type":"linear","range":"width","zero":false,"domain":{"data":"bf701b1c-6afd-4be4-92ab-b9d9ec32a288","field":"data.x"}},{"name":"y","type":"linear","range":"height","nice":true,"zero":false,"domain":{"data":"bf701b1c-6afd-4be4-92ab-b9d9ec32a288","field":"data.y"}}],"marks":[{"type":"line","from":{"data":"bf701b1c-6afd-4be4-92ab-b9d9ec32a288"},"properties":{"enter":{"x":{"scale":"x","field":"data.x"},"y":{"scale":"y","field":"data.y"},"stroke":{"value":"#FF29D2"},"strokeWidth":{"value":2},"strokeOpacity":{"value":1}}}}],"data":[{"name":"bf701b1c-6afd-4be4-92ab-b9d9ec32a288","values":[{"x":0,"y":4.09},{"x":1,"y":3.395},{"x":2,"y":3.715},{"x":3,"y":3.215},{"x":4,"y":3.32},{"x":5,"y":2.79},{"x":6,"y":3.565},{"x":7,"y":2.825},{"x":8,"y":3.64},{"x":9,"y":3.77},{"x":10,"y":3.19},{"x":11,"y":2.08},{"x":12,"y":1.56},{"x":13,"y":1.565},{"x":14,"y":1.59},{"x":15,"y":1.785},{"x":16,"y":1.84},{"x":17,"y":1.955},{"x":18,"y":2.0},{"x":19,"y":1.935},{"x":20,"y":2.375},{"x":21,"y":2.655},{"x":22,"y":2.865},{"x":23,"y":2.81},{"x":24,"y":2.685},{"x":25,"y":2.535},{"x":26,"y":2.79},{"x":27,"y":2.845},{"x":28,"y":2.81},{"x":29,"y":2.795},{"x":30,"y":2.685},{"x":31,"y":2.675},{"x":32,"y":2.855},{"x":33,"y":2.66},{"x":34,"y":2.52},{"x":35,"y":2.76},{"x":36,"y":2.755},{"x":37,"y":2.6},{"x":38,"y":2.73},{"x":39,"y":2.795},{"x":40,"y":2.725},{"x":41,"y":3.165},{"x":42,"y":2.945},{"x":43,"y":2.785},{"x":44,"y":3.05},{"x":45,"y":2.87},{"x":46,"y":3.18},{"x":47,"y":3.12},{"x":48,"y":3.07},{"x":49,"y":3.365},{"x":50,"y":2.935},{"x":51,"y":2.625},{"x":52,"y":2.63},{"x":53,"y":2.755},{"x":54,"y":2.74},{"x":55,"y":2.725},{"x":56,"y":3.16},{"x":57,"y":2.975},{"x":58,"y":2.575},{"x":59,"y":3.015},{"x":60,"y":2.93},{"x":61,"y":2.765},{"x":62,"y":2.715},{"x":63,"y":2.99},{"x":64,"y":2.98},{"x":65,"y":3.07},{"x":66,"y":3.19},{"x":67,"y":3.195},{"x":68,"y":3.075},{"x":69,"y":3.23},{"x":70,"y":3.085},{"x":71,"y":3.085},{"x":72,"y":3.14},{"x":73,"y":3.495},{"x":74,"y":3.15},{"x":75,"y":3.29},{"x":76,"y":3.405},{"x":77,"y":3.145},{"x":78,"y":3.185},{"x":79,"y":3.33},{"x":80,"y":3.18},{"x":81,"y":2.85},{"x":82,"y":2.975},{"x":83,"y":2.69},{"x":84,"y":3.18},{"x":85,"y":3.375},{"x":86,"y":2.475},{"x":87,"y":2.405},{"x":88,"y":2.59},{"x":89,"y":2.625},{"x":90,"y":2.615},{"x":91,"y":2.66},{"x":92,"y":2.605},{"x":93,"y":3.06},{"x":94,"y":2.82},{"x":95,"y":2.845},{"x":96,"y":2.885},{"x":97,"y":2.785},{"x":98,"y":2.865},{"x":99,"y":3.155},{"x":100,"y":2.84},{"x":101,"y":2.985},{"x":102,"y":3.085},{"x":103,"y":2.375},{"x":104,"y":2.775},{"x":105,"y":2.47},{"x":106,"y":2.77},{"x":107,"y":2.76},{"x":108,"y":2.855},{"x":109,"y":3.065},{"x":110,"y":3.1},{"x":111,"y":2.645},{"x":112,"y":2.41},{"x":113,"y":2.565},{"x":114,"y":2.695},{"x":115,"y":2.6},{"x":116,"y":2.32},{"x":117,"y":2.355},{"x":118,"y":2.085},{"x":119,"y":1.985},{"x":120,"y":2.05},{"x":121,"y":1.955},{"x":122,"y":1.91},{"x":123,"y":1.815},{"x":124,"y":2.025},{"x":125,"y":1.97},{"x":126,"y":2.11},{"x":127,"y":2.435},{"x":128,"y":2.58},{"x":129,"y":2.425},{"x":130,"y":2.4},{"x":131,"y":2.63},{"x":132,"y":2.67},{"x":133,"y":2.63},{"x":134,"y":2.685},{"x":135,"y":2.735},{"x":136,"y":2.66},{"x":137,"y":2.785},{"x":138,"y":2.87},{"x":139,"y":2.345},{"x":140,"y":2.435},{"x":141,"y":2.405},{"x":142,"y":2.645},{"x":143,"y":2.855},{"x":144,"y":2.61},{"x":145,"y":2.775},{"x":146,"y":2.81},{"x":147,"y":2.805},{"x":148,"y":3.16},{"x":149,"y":3.18},{"x":150,"y":3.075},{"x":151,"y":2.915},{"x":152,"y":2.25},{"x":153,"y":2.56},{"x":154,"y":2.795},{"x":155,"y":2.835},{"x":156,"y":2.915},{"x":157,"y":2.795},{"x":158,"y":2.785},{"x":159,"y":2.785},{"x":160,"y":2.625},{"x":161,"y":3.095},{"x":162,"y":2.565},{"x":163,"y":2.655},{"x":164,"y":2.875},{"x":165,"y":2.985},{"x":166,"y":2.53},{"x":167,"y":2.52},{"x":168,"y":2.785},{"x":169,"y":3.005},{"x":170,"y":3.06},{"x":171,"y":3.265},{"x":172,"y":3.16},{"x":173,"y":3.135},{"x":174,"y":3.225},{"x":175,"y":3.125},{"x":176,"y":3.005},{"x":177,"y":3.095},{"x":178,"y":3.41},{"x":179,"y":3.12},{"x":180,"y":3.44},{"x":181,"y":3.0},{"x":182,"y":3.17},{"x":183,"y":3.335},{"x":184,"y":3.27},{"x":185,"y":3.48},{"x":186,"y":3.41},{"x":187,"y":3.425},{"x":188,"y":3.65},{"x":189,"y":3.11},{"x":190,"y":3.42},{"x":191,"y":3.465},{"x":192,"y":3.32},{"x":193,"y":3.4},{"x":194,"y":3.38},{"x":195,"y":3.125},{"x":196,"y":3.69},{"x":197,"y":3.455},{"x":198,"y":3.25},{"x":199,"y":3.285}]}],"width":400,"height":247.2187957763672,"padding":{"bottom":20,"top":10,"right":10,"left":50}},"value":"#gorilla_repl.vega.VegaView{:content {:axes [{:scale \"x\", :type \"x\"} {:scale \"y\", :type \"y\"}], :scales [{:name \"x\", :type \"linear\", :range \"width\", :zero false, :domain {:data \"bf701b1c-6afd-4be4-92ab-b9d9ec32a288\", :field \"data.x\"}} {:name \"y\", :type \"linear\", :range \"height\", :nice true, :zero false, :domain {:data \"bf701b1c-6afd-4be4-92ab-b9d9ec32a288\", :field \"data.y\"}}], :marks [{:type \"line\", :from {:data \"bf701b1c-6afd-4be4-92ab-b9d9ec32a288\"}, :properties {:enter {:x {:scale \"x\", :field \"data.x\"}, :y {:scale \"y\", :field \"data.y\"}, :stroke {:value \"#FF29D2\"}, :strokeWidth {:value 2}, :strokeOpacity {:value 1}}}}], :data [{:name \"bf701b1c-6afd-4be4-92ab-b9d9ec32a288\", :values ({:x 0, :y 4.09} {:x 1, :y 3.395} {:x 2, :y 3.715} {:x 3, :y 3.215} {:x 4, :y 3.32} {:x 5, :y 2.79} {:x 6, :y 3.565} {:x 7, :y 2.825} {:x 8, :y 3.64} {:x 9, :y 3.77} {:x 10, :y 3.19} {:x 11, :y 2.08} {:x 12, :y 1.56} {:x 13, :y 1.565} {:x 14, :y 1.59} {:x 15, :y 1.785} {:x 16, :y 1.84} {:x 17, :y 1.955} {:x 18, :y 2.0} {:x 19, :y 1.935} {:x 20, :y 2.375} {:x 21, :y 2.655} {:x 22, :y 2.865} {:x 23, :y 2.81} {:x 24, :y 2.685} {:x 25, :y 2.535} {:x 26, :y 2.79} {:x 27, :y 2.845} {:x 28, :y 2.81} {:x 29, :y 2.795} {:x 30, :y 2.685} {:x 31, :y 2.675} {:x 32, :y 2.855} {:x 33, :y 2.66} {:x 34, :y 2.52} {:x 35, :y 2.76} {:x 36, :y 2.755} {:x 37, :y 2.6} {:x 38, :y 2.73} {:x 39, :y 2.795} {:x 40, :y 2.725} {:x 41, :y 3.165} {:x 42, :y 2.945} {:x 43, :y 2.785} {:x 44, :y 3.05} {:x 45, :y 2.87} {:x 46, :y 3.18} {:x 47, :y 3.12} {:x 48, :y 3.07} {:x 49, :y 3.365} {:x 50, :y 2.935} {:x 51, :y 2.625} {:x 52, :y 2.63} {:x 53, :y 2.755} {:x 54, :y 2.74} {:x 55, :y 2.725} {:x 56, :y 3.16} {:x 57, :y 2.975} {:x 58, :y 2.575} {:x 59, :y 3.015} {:x 60, :y 2.93} {:x 61, :y 2.765} {:x 62, :y 2.715} {:x 63, :y 2.99} {:x 64, :y 2.98} {:x 65, :y 3.07} {:x 66, :y 3.19} {:x 67, :y 3.195} {:x 68, :y 3.075} {:x 69, :y 3.23} {:x 70, :y 3.085} {:x 71, :y 3.085} {:x 72, :y 3.14} {:x 73, :y 3.495} {:x 74, :y 3.15} {:x 75, :y 3.29} {:x 76, :y 3.405} {:x 77, :y 3.145} {:x 78, :y 3.185} {:x 79, :y 3.33} {:x 80, :y 3.18} {:x 81, :y 2.85} {:x 82, :y 2.975} {:x 83, :y 2.69} {:x 84, :y 3.18} {:x 85, :y 3.375} {:x 86, :y 2.475} {:x 87, :y 2.405} {:x 88, :y 2.59} {:x 89, :y 2.625} {:x 90, :y 2.615} {:x 91, :y 2.66} {:x 92, :y 2.605} {:x 93, :y 3.06} {:x 94, :y 2.82} {:x 95, :y 2.845} {:x 96, :y 2.885} {:x 97, :y 2.785} {:x 98, :y 2.865} {:x 99, :y 3.155} {:x 100, :y 2.84} {:x 101, :y 2.985} {:x 102, :y 3.085} {:x 103, :y 2.375} {:x 104, :y 2.775} {:x 105, :y 2.47} {:x 106, :y 2.77} {:x 107, :y 2.76} {:x 108, :y 2.855} {:x 109, :y 3.065} {:x 110, :y 3.1} {:x 111, :y 2.645} {:x 112, :y 2.41} {:x 113, :y 2.565} {:x 114, :y 2.695} {:x 115, :y 2.6} {:x 116, :y 2.32} {:x 117, :y 2.355} {:x 118, :y 2.085} {:x 119, :y 1.985} {:x 120, :y 2.05} {:x 121, :y 1.955} {:x 122, :y 1.91} {:x 123, :y 1.815} {:x 124, :y 2.025} {:x 125, :y 1.97} {:x 126, :y 2.11} {:x 127, :y 2.435} {:x 128, :y 2.58} {:x 129, :y 2.425} {:x 130, :y 2.4} {:x 131, :y 2.63} {:x 132, :y 2.67} {:x 133, :y 2.63} {:x 134, :y 2.685} {:x 135, :y 2.735} {:x 136, :y 2.66} {:x 137, :y 2.785} {:x 138, :y 2.87} {:x 139, :y 2.345} {:x 140, :y 2.435} {:x 141, :y 2.405} {:x 142, :y 2.645} {:x 143, :y 2.855} {:x 144, :y 2.61} {:x 145, :y 2.775} {:x 146, :y 2.81} {:x 147, :y 2.805} {:x 148, :y 3.16} {:x 149, :y 3.18} {:x 150, :y 3.075} {:x 151, :y 2.915} {:x 152, :y 2.25} {:x 153, :y 2.56} {:x 154, :y 2.795} {:x 155, :y 2.835} {:x 156, :y 2.915} {:x 157, :y 2.795} {:x 158, :y 2.785} {:x 159, :y 2.785} {:x 160, :y 2.625} {:x 161, :y 3.095} {:x 162, :y 2.565} {:x 163, :y 2.655} {:x 164, :y 2.875} {:x 165, :y 2.985} {:x 166, :y 2.53} {:x 167, :y 2.52} {:x 168, :y 2.785} {:x 169, :y 3.005} {:x 170, :y 3.06} {:x 171, :y 3.265} {:x 172, :y 3.16} {:x 173, :y 3.135} {:x 174, :y 3.225} {:x 175, :y 3.125} {:x 176, :y 3.005} {:x 177, :y 3.095} {:x 178, :y 3.41} {:x 179, :y 3.12} {:x 180, :y 3.44} {:x 181, :y 3.0} {:x 182, :y 3.17} {:x 183, :y 3.335} {:x 184, :y 3.27} {:x 185, :y 3.48} {:x 186, :y 3.41} {:x 187, :y 3.425} {:x 188, :y 3.65} {:x 189, :y 3.11} {:x 190, :y 3.42} {:x 191, :y 3.465} {:x 192, :y 3.32} {:x 193, :y 3.4} {:x 194, :y 3.38} {:x 195, :y 3.125} {:x 196, :y 3.69} {:x 197, :y 3.455} {:x 198, :y 3.25} {:x 199, :y 3.285})}], :width 400, :height 247.2188, :padding {:bottom 20, :top 10, :right 10, :left 50}}}"}
;; <=

;; @@

;; @@

;; @@

;; @@
