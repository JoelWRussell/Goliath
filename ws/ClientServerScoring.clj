;; gorilla-repl.fileformat = 1

;; **
;;; # ClientServerScoring
;;; This worksheet demonstrates networked scoring and the finding of coefficients. The LagrangeServer needs to be running somewhere - maybe on a different machine. There also needs to be at least one WorkerClient also possibly on a different machine. Or several on different machines, but they must have mma on them and have JLink etc. This worksheet will first of all contact a www address to find the location of the server - which saves remembering ip:port, it will make contact with the server and then control the group of worker computers via commands like InitMathKernel, NewData, NewZeitgeist, GetScores, GetCoefficients. The LagrangeServer will sort out sharing out the task etc and will just return the results. 
;; **

;; **
;;; # GorillaClient
;;; This is sends data and results from this worksheet to the LagrangeServer - it knows how to talk to this worksheet and also how to find LagrangeServer on the web and send data etc.
;; **

;; **
;;; #LagrangeServer
;;; This manages all of the WorkerClients which have mma installed. Slices up the task of scoring and collects results.
;; **

;; **
;;; #WorkerClient
;;; This has mma and jlink etc and actually does the calculation. It knows how to find LagrangeServer on the web and how to talk to it.
;; **

;; @@
(ns combative-canopy
  (:require [gorilla-plot.core :as plot]))
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@
(import '[com.lagrangianmining.GorillaClient])
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@
(def client (com.lagrangianmining.GorillaClient.)) ;;GorillaClient finds LagrangeServer on the web and shuttles commands to it and gets data from it.
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;combative-canopy/client</span>","value":"#'combative-canopy/client"}
;; <=

;; @@
(.InitMathKernel client);;Create a MMA Kernel on all of the WorkerClients which are registered with LagrangeServer.
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>true</span>","value":"true"}
;; <=

;; @@
(.NewData client "resources/sho_coupled_0.1_sim.csv") ;;transfers the experimental data across to all of the WorkerClients
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>true</span>","value":"true"}
;; <=

;; @@
(.PrepareData client true 2 0.1) ;;tells the WorkerClients to load the PrepareDataCompact.m script file
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>true</span>","value":"true"}
;; <=

;; @@
;;(.ResetLagrangeServer client) ;;this will destroy and restart the remote LagrangeServer
;;all of the workers will reconnect
;; @@

;; @@
(def polyLength (into-array Integer/TYPE [20]))
(def poly (into-array Integer/TYPE [2 0 0 0 0 2 0 0 1 1 0 0 0 0 2 0 0 0 0 2]))
(def deltat 0.1)
(def df 2)

;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;combative-canopy/df</span>","value":"#'combative-canopy/df"}
;; <=

;; @@
(.NewZeitgeist client 1 polyLength poly df deltat );;score the zeitgeist with the netwrok
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>true</span>","value":"true"}
;; <=

;; @@
(into [] (.GetScores client))
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>-22.806634470647438</span>","value":"-22.806634470647438"}],"value":"[-22.806634470647438]"}
;; <=

;; @@
(into [] (.GetCoefficients client))
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>0.18854088335408387</span>","value":"0.18854088335408387"},{"type":"html","content":"<span class='clj-double'>0.09427039337015511</span>","value":"0.09427039337015511"},{"type":"html","content":"<span class='clj-double'>-0.188540811738729</span>","value":"-0.188540811738729"},{"type":"html","content":"<span class='clj-double'>-0.09427045853381483</span>","value":"-0.09427045853381483"},{"type":"html","content":"<span class='clj-double'>-0.09427042790594096</span>","value":"-0.09427042790594096"}],"value":"[0.18854088335408387 0.09427039337015511 -0.188540811738729 -0.09427045853381483 -0.09427042790594096]"}
;; <=

;; @@
(.DestroyMathKernel client) ;;closes the mathkernel on all of the client computers
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>true</span>","value":"true"}
;; <=

;; @@
(.DisconnectServer client);;breaks the connection with LagrangeServer - another computer could make contact at this point
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-nil'>nil</span>","value":"nil"}
;; <=

;; @@

;; @@
