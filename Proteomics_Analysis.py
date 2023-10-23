https://github.com/Nesvilab/MSFragger

Building from scratch
Update build version:
The version of the build is stored in 3 separate places:

File: MSFragger-GUI/src/umich/msfragger/gui/Bundle.properties
Property: msfragger.gui.version
File: MSFragger-GUI/build.gradle
Property: version
File: MSFragger-GUI/src/umich/msfragger/gui/Bundle.properties 
Property: msfragger.gui.version
Build:
You don't need to have Gradle installed, the Gradle wrapper included in this repository will be used. From the root directory of the repository issue the following commands:

cd ./MSFragger-GUI
./gradlew makeReleaseZipNoJre
or use this version to build with Java Runtime (for Windows only):

cd ./MSFragger-GUI
./gradlew makeReleaseZipWithJre
The .zip output will be in MSFragger-GUI/build/github-release.
