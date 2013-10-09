package org.clothocad.core.testers;

import java.net.InetAddress;

public class PrintLocalhost {
  public static void main(String[] args) throws Exception{
    InetAddress localhost = InetAddress.getLocalHost();
    String host = localhost.getHostAddress();
    System.out.println("host: "+host);
  }
}