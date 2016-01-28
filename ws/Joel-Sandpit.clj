;; gorilla-repl.fileformat = 1

;; **
;;; # Joel-Sandpit
;;; 
;;; A place for trying out the cross-over and mutate functions.
;;; 
;;; The score function can be called with ts/poly-score
;;; 
;; **

;; @@
(ns goliath
  (:require [gorilla-plot.core :as plot]
            [score.regression-test-score :as ts]))
;; @@

;; @@
(ts/poly-score [1 3])
;; @@

;; **
;;; The cross over function must take two polynomials and produce two child polynomials. 
;; **

;; @@
; Random number generator (excluding 0)

(defn random-num
  []
  "Generates a random number excluding the value 0"
  (let [x (rand)]
    (if (= x 0)
       (+ x (rand))) 
      x))


;need to add a final check if x is 0.
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;goliath/random-num</span>","value":"#'goliath/random-num"}
;; <=

;; @@
;Take to vectors    [[0 2]] [[2 3] [2 3]]]    [[0 2] [1 3] [3 5]]
; Created two empty chidren of length equal to partents.
; Assign genes from parent one to child one with probability X
; Else assign gene from parent one to child to poraility 1-X
; So long as the child genes do not exceed thier parental legnths add genes from parent 1 to offspring with the same probaiblity of.



(defn splice 
  [parent-gene child]
  (conj child parent-gene))
g
;if in both put it into both children 
;if length of child is the same or greater than the lenght of parent put into other 



(defn cross-over [poly1 poly2]
  (let [c1 #{}
    	c2 #{}
        rnd (rand)] 
    
    (if (< (rand) rnd)
    
(mapv (fn [poly1 poly2] (conj [(vec x)] (vec y))) poly1 poly2)
   ))

(contains? [[1 3] [1 2]] [1 3])


; the first genotype from parent 1 goes to child 1 with probability x
; else genotype one goes to child 2.
; offspring should have the same length as their parents.
; distribute polynimal 1 first



;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>false</span>","value":"false"}
;; <=

;; @@
(def x [[1 2] [3 5] [7 2] [6 0]] )
(def y [[3 4] [2 4] [3 5] [1 5]] )

(def c1 [])
(def c2 [])     

;(mapv (fn [p1 p2 c1 c2]
 ;       (if (<= rnd (rand))
  ;        bn        
   ;       )
  ;(defn child-allocation
  ;[parent child-1 child-2]  
  ;)

        
(defn add-gene 
   [parent i child-1 child-2]
  "Adds the ith gene of parent to child"
  (println parent)
  (println child-1)
  (println child-2)
  (if (contains? (set child-1) (parent i))
  (concat child-1 (parent i))
  (concat child-2 (parent i)))
   parent child-1 child-2)

(add-gene x 1 c1 c2)



(if (contains? (set child-1) (x 1) ))


  

(defn add-gene 
   [parent i child-1 child-2]
  "Adds the ith gene of parent to child"
  (if (contains? (set child-1) (parent i))
  (concat child-1 (parent i))
  (concat child-2 (parent i)))
   parent child-1 child-2)
;; @@
;; ->
;;; [[1 2] [3 5] [7 2] [6 0]]
;;; []
;;; []
;;; 
;; <-
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/add-gene</span>","value":"#'user/add-gene"}
;; <=

;; @@
;(def x [[1 2] [2 4] [4 1]])
;(def child-x [])

(defn add-gene 
   [parent i child]
  "adds allele i from partent to child"
  (conj child (parent i)))

;(def child-x (add-gene x 1 child-x))

;child-x
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-long'>4</span>","value":"4"}],"value":"[2 4]"}],"value":"[[2 4]]"}
;; <=

;; @@
(defn conditional-add 
  [random-number p1 c1 c2 i]
  (if (<= random-number (rand))
       (add-gene p1 i c1)
       (add-gene p1 i c2)))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/conditional-add</span>","value":"#'user/conditional-add"}
;; <=

;; @@
(def random_num 0.5)
(def child1 [])
(def child2 [])
(def parent [[1 0] [1 1] [3 4]])
  
(mapv (fn [x] (if (<= random_num (rand))
          "child1"
          "child2")) parent)

(reduce (fn [x y] 
          (if (<= y (rand))
          (conj child1 x)
          (conj child2 x)) random_num parent)
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-string'>&quot;child1&quot;</span>","value":"\"child1\""},{"type":"html","content":"<span class='clj-string'>&quot;child2&quot;</span>","value":"\"child2\""},{"type":"html","content":"<span class='clj-string'>&quot;child2&quot;</span>","value":"\"child2\""}],"value":"[\"child1\" \"child2\" \"child2\"]"}
;; <=

;; @@
(if (contains? parent2 x))
	 () 
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-long'>25</span>","value":"25"}
;; <=

;; @@
(defn distribute-parent-1 
  [parent-1 child-1 child-2]
  "Distributes the alleles of parent 1 to children with random probability "
  (let [rnd-num (rand) 
        p1 parent-1 
        c1 child-1
        c2 child-2
        length_p (count p1)
        value 0.0]
    
   
    (if (< length_p (+ value 1))
      [c1 c2]  
      (recur (conditional-add rnd-num p1 c1 c2 value) (inc value)))
    [p1 c1 c2])
  
  
 (func)


 (mapv (if )
;; @@

;; @@
(def len 4)
(def value 0)

(if (< len (+ 0 value))
       value
    (recur (inc value)))   
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/value</span>","value":"#'user/value"}
;; <=

;; @@
(mapv (fn [x y] (conj [(vec x)] (vec y))) x y)

;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-long'>3</span>","value":"3"}],"value":"[1 3]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>2</span>","value":"2"},{"type":"html","content":"<span class='clj-long'>4</span>","value":"4"}],"value":"[2 4]"}],"value":"[[1 3] [2 4]]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>3</span>","value":"3"},{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"}],"value":"[3 1]"},{"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-long'>1</span>","value":"1"},{"type":"html","content":"<span class='clj-long'>0</span>","value":"0"}],"value":"[1 0]"}],"value":"[[3 1] [1 0]]"}],"value":"[[[1 3] [2 4]] [[3 1] [1 0]]]"}
;; <=

;; @@

;; @@

;; @@
(def parent1 [[1 2] [0 1] [2 0]])
(def parent2 [[1 4] [1 0] [3 0]])

(def offspring {:child1 []  :child2 []})

;(conj (offspring :child1) [1 0])


(defn joel
  [p1 p2]
  (loop []))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;user/offspring</span>","value":"#'user/offspring"}
;; <=

;; @@
(conj offspring)
;; @@
