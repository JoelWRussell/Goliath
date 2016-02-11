# Goliath

A hybrid polynomial and genetic algorithm search for Lagrangians. 


## Installing JLink

Before running Goliath for the first time you will need to install the JLink jar, which allows Java code to call the Mathematica kernel, into your local maven repository. This will allow Leiningen to find the library. To do this, first find out where the JLink.jar file lives on your computer. On a Mac it's likely to be somewhere like `/Applications/Mathematica.app/Contents/SystemFiles/Links/JLink/JLink.jar`. Then, from within the Goliath directory, run the command:
```
lein localrepo install /path/to/JLink.jar jlink/jlink 1.0.0
```
This will install the jar into the Maven repo, so that it can be referenced frome the `project.clj` file as a dependency (`[jlink "1.0.0"]). Note that this dependency is already there - you don't need to add it - you just need to run the above command so Leiningen knows how to resolve it.

The next thing you need to do is to make sure that JLink can find the native library that talks to the kernel. The way to do this is nasty, forced upon us by the fact that Wolfram really don't like to make it easy to use anything other than Mathematica. The way to do this is to copy the native libs to the root of the project, which is on the java.library.path by default. On a Mac you'll want to do something like this, from within the Goliath directory:
```
cp /Applications/Mathematica.app/Contents/SystemFiles/Links/JLink/SystemFiles/Libraries/MacOSX-x86-64/* .
```
Windows will be similar, but the path will be different. Don't commit these library files to git. There are `.gitignore` entries that should stop them, but keep an eye out anyway.
