FROM maven:3.3.3-jdk-8
MAINTAINER David Tran <djtran@bu.edu>

COPY . clotho/
WORKDIR clotho/


CMD ["mvn", "'-Dexec.args=-classpath %classpath org.clothocad.core.util.ClothoTestEnvironment' -Dexec.executable=/usr/bin/java -Dexec.classpathScope=test org.codehaus.mojo:exec-maven-plugin:1.2.1:exec"]

