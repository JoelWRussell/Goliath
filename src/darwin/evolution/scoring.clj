;
; This file is part of darwin.
;
; Copyright (C) 2014-, Imperial College, London, All rights reserved.
;
; Contributors: Jony Hudson
;
; Released under the MIT license..
;

(ns darwin.evolution.scoring

)


;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;this function was originally called to score each one individually
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn update-individual-scores
  "Update the scores for an individual. The score functions are given as a map of functions:
  each function will be applied to the individual's genotype and its result stored on the individual,
  under the function's key."
  [score-funcs individual]
  (merge individual
         (into {} (map (fn [s] [(first s) ((second s) (:genotype individual))]) score-funcs))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;this function scores a the genome of the whole population via the network
;;assume that all of the appropriate network connections have been set up
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn update-scores
  "Update the scores for each individual in the given list. See above for how the score functions are
   specified."
  [individuals score-funcs]
 (let [gen (map :genotype individuals)
	res (user/networkScore gen)
	xx (map (fn [x1 x2] (hash-map :complexity (count (flatten x1)) :error x2 )) gen res)
	yy (map #(merge %1 %2) individuals xx)]
	yy
))
