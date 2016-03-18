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
(.NewData client "resources/mma_double.csv") ;;transfers the experimental data across to all of the WorkerClients
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>true</span>","value":"true"}
;; <=

;; @@
(.PrepareData client false 2 0.1) ;;tells the WorkerClients to load the PrepareDataCompact.m script file
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>true</span>","value":"true"}
;; <=

;; @@
(def polyLength (into-array Integer/TYPE [4 4 4]))
(def poly (into-array Integer/TYPE [1 1 1 1  2 2 2 2 1 2 1 2]))
(def deltat 0.1)
(def df 2)

;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-var'>#&#x27;combative-canopy/df</span>","value":"#'combative-canopy/df"}
;; <=

;; @@
(.NewZeitgeist client 3 polyLength poly df deltat );;score the zeitgeist with the netwrok
;; @@
;; =>
;;; {"type":"html","content":"<span class='clj-unkown'>true</span>","value":"true"}
;; <=

;; @@
(into [] (.GetScores client))
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>14.072445026423615</span>","value":"14.072445026423615"},{"type":"html","content":"<span class='clj-double'>22.928462152877966</span>","value":"22.928462152877966"},{"type":"html","content":"<span class='clj-double'>18.157463275646037</span>","value":"18.157463275646037"}],"value":"[14.072445026423615 22.928462152877966 18.157463275646037]"}
;; <=

;; @@
(into [] (.GetCoefficients client))
;; @@
;; =>
;;; {"type":"list-like","open":"<span class='clj-vector'>[</span>","close":"<span class='clj-vector'>]</span>","separator":" ","items":[{"type":"html","content":"<span class='clj-double'>8.129417044365265</span>","value":"8.129417044365265"},{"type":"html","content":"<span class='clj-double'>261.3891625397806</span>","value":"261.3891625397806"},{"type":"html","content":"<span class='clj-double'>39.27242908652074</span>","value":"39.27242908652074"}],"value":"[8.129417044365265 261.3891625397806 39.27242908652074]"}
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
