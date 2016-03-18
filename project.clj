;;;; This file is part of flow. Copyright (C) 2014-, Jony Hudson.
;;;;
;;;; Released under the MIT license.
;;;;

(defproject Goliath "1.0.0"
  :url "https://github.com/JoelWRussell/Goliath"
  :license {:name "MIT"}
  :dependencies [[org.clojure/clojure "1.7.0-alpha5"]
                 [gorilla-plot "0.1.3"]
                 [algebolic "1.0.0"]
                 [darwin "1.0.0"]
                 [semantic-csv "0.1.0-alpha1"]
                 [criterium "0.4.3"]
                 [incanter-gorilla "0.1.0"]
                 [jlink "1.0.0"]
		 [commons-codec-1.9 "1.0.0"]
		 [commons-logging-1.2 "1.0.0"]
		 [fluent-hc-4.5.2 "1.0.0"]
		 [httpclient-4.5.2 "1.0.0"]
		 [httpclient-cache-4.5.2 "1.0.0"]
		 [httpcore-4.4.4 "1.0.0"]
		 [httpmime-4.5.2 "1.0.0"]
		 [jna-4.1.0 "1.0.0"]
		 [jna-platform-4.1.0 "1.0.0"]

]
		
  :java-source-paths ["java"]
  :plugins [[lein-gorilla "0.3.5"]
            [lein-localrepo "0.5.3"]]
  :jvm-opts ^:replace ["-server"
                       ;;"-XX:+AggressiveOpts"
                       ;;"-XX:+UseFastAccessorMethods"
                       ;;"-XX:+UseCompressedOops"
                       "-Xmx1g"]
)
