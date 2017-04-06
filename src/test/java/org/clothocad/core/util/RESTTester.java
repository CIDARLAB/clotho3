/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.MalformedURLException;
import java.net.ProtocolException;
import java.net.URL;
import java.security.KeyManagementException;
import java.security.NoSuchAlgorithmException;
import java.util.HashMap;
import java.util.Map;
import javax.net.ssl.HttpsURLConnection;
import javax.net.ssl.SSLContext;
import javax.net.ssl.TrustManager;
import javax.net.ssl.X509TrustManager;
import org.junit.Test;

/**
 *
 * @author David
 */
public class RESTTester {

    private String url = "https://localhost:8443/data/post";

    TrustManager[] trustAllCerts = new TrustManager[]{
        new X509TrustManager() {
            public java.security.cert.X509Certificate[] getAcceptedIssuers() {
                return null;
            }

            public void checkClientTrusted(
                    java.security.cert.X509Certificate[] certs, String authType) {
            }

            public void checkServerTrusted(
                    java.security.cert.X509Certificate[] certs, String authType) {
            }
        }
    };

    public String HTTPReq(URL url, String jsonString, String verb) throws ProtocolException, IOException, KeyManagementException, NoSuchAlgorithmException {

        SSLContext sc = SSLContext.getInstance("SSL");
        sc.init(null, trustAllCerts, new java.security.SecureRandom());
        HttpsURLConnection.setDefaultSSLSocketFactory(sc.getSocketFactory());

        HttpsURLConnection conn = (HttpsURLConnection) url.openConnection();
        conn.setDoInput(true);
        conn.setRequestMethod(verb);
        conn.setRequestProperty("Content-Type", "application/json");
        conn.setRequestProperty("Cache-Control", "no-cache");
        conn.setInstanceFollowRedirects(false);
        if (!verb.equals("GET")) {
            conn.setDoOutput(true);
            OutputStream os = conn.getOutputStream();
            os.write(jsonString.getBytes());
            os.flush();
        }

        int responseCode = conn.getResponseCode();

        if (responseCode != 400 && responseCode != 404 && responseCode != 500) {
//            System.out.println("SUCCESS!");

            //print result
            BufferedReader br = new BufferedReader(new InputStreamReader((conn.getInputStream())));

            String output;
            String alloutput = "";
            while ((output = br.readLine()) != null) {
                alloutput += output;
            }
            conn.disconnect();
            return alloutput;
        } else {
            conn.disconnect();
            return "ERROR " + responseCode;
        }

    }

//    @Test
    public void testCreateUser() throws MalformedURLException, IOException, KeyManagementException, NoSuchAlgorithmException {
        System.out.println("Testing Create User");
        String jsonString = "{'username':'jsmith1','password':'asdf'}";
        URL url = new URL(this.url + "/create/user");

        String output = HTTPReq(url, jsonString, "POST");

        System.out.println(output);
    }
    
    @Test
    public void testCreateSequence() throws MalformedURLException, IOException, KeyManagementException, NoSuchAlgorithmException {
        System.out.println("Testing Create Sequence");
        String jsonString = "{'username':'jsmith','objectName':'Test Sequence','sequence':'ata'}";
        URL url = new URL(this.url + "/create/sequence/");

        String output = HTTPReq(url, jsonString, "POST");

        System.out.println(output);
    }

    @Test
    public void testCreatePart() throws MalformedURLException, IOException, KeyManagementException, NoSuchAlgorithmException {

        String jsonString = "{'username':'jsmith','objectName':'Test Sequence','sequence':'ata'}";
        URL url = new URL(this.url + "/create/sequence");
        String seqId = HTTPReq(url, jsonString, "POST");

        System.out.println("Testing Create Part");
        jsonString = "{'username':'jsmith','objectName':'Test Part', 'id':'" + seqId + "'}";
        url = new URL(this.url + "/create/part");

        String output = HTTPReq(url, jsonString, "POST");

        System.out.println(output);
    }

    @Test
    public void testGetByName() throws MalformedURLException, IOException, KeyManagementException, NoSuchAlgorithmException {
        System.out.println("Testing Get Sequence by Name");

        String jsonString = "{'username':'jsmith','objectName':'TestSequence','sequence':'ata'}";
        URL url = new URL(this.url + "/create/sequence");
        String seqId = HTTPReq(url, jsonString, "POST");

        url = new URL("https://localhost:8443/data/get/getByName/TestSequence/");

        String output = HTTPReq(url, "", "GET");

        System.out.println(output);
    }

    @Test
    public void testGetById() throws MalformedURLException, IOException, KeyManagementException, NoSuchAlgorithmException {
        System.out.println("Testing Get By Id");

        String jsonString = "{'username':'jsmith','objectName':'Test Sequence','sequence':'ata'}";
        URL url = new URL(this.url + "/create/sequence");
        String seqId = HTTPReq(url, jsonString, "POST");

        url = new URL("https://localhost:8443/data/get/getById/" + seqId);

        String output = HTTPReq(url, "", "GET");

        System.out.println(output);

    }

    @Test
    public void testSet() throws MalformedURLException, IOException, KeyManagementException, NoSuchAlgorithmException {

<<<<<<< HEAD
        String jsonString = "{'name':'jsmith','objectName':'Test Sequence','sequence':'ata'}";
=======
        String jsonString = "{'username':'jsmith','objectName':'Test Sequence','sequence':'ata'}";
>>>>>>> d24bd7443a2e06fcae28c1dc5243f18cdc592221
        URL url = new URL(this.url + "/create/sequence");
        String seqId = HTTPReq(url, jsonString, "POST");

        System.out.println("Testing Set");
<<<<<<< HEAD
        jsonString = "{'name':'jsmith', 'sequence' : 'atatatatatatat','id' : '" + seqId + "a'}";
=======
        jsonString = "{'username':'jsmith','sequence' : 'atatatatatatat','id' : '" + seqId + "'}";
>>>>>>> d24bd7443a2e06fcae28c1dc5243f18cdc592221
        url = new URL(this.url + "/set");

        String output = HTTPReq(url, jsonString, "PUT");

        System.out.println(output);
    }

    @Test
    public void testDelete() throws MalformedURLException, IOException, KeyManagementException, NoSuchAlgorithmException {
        System.out.println("Testing Delete Sequence");

        String jsonString = "{'username':'jsmith','objectName':'Test Sequence','sequence':'ata'}";
        URL url = new URL(this.url + "/create/sequence");

        String sequenceId = HTTPReq(url, jsonString, "POST");

        jsonString = "{'username':'jsmith','id':'" + sequenceId + "'}";
        url = new URL(this.url + "/delete/delete");

        String output = HTTPReq(url, jsonString, "DELETE");

        System.out.println(output);
    }

    @Test
    public void testConveniencePart() throws MalformedURLException, IOException, KeyManagementException, NoSuchAlgorithmException {
        System.out.println("Testing Create Convenience Part");

        String jsonString = "{'username':'jsmith','objectName':'Test Convenience Part','sequence':'tccctatcagtgatagagattgacatccctatcagtgatagagatactgagcac', 'role':'GENE'}";
        URL url = new URL(this.url + "/create/conveniencePart/");

        String output = HTTPReq(url, jsonString, "POST");

        System.out.println(output);
    }
    
    @Test
    public void testConvenienceDevice() throws MalformedURLException, IOException, KeyManagementException, NoSuchAlgorithmException {
        System.out.println("Testing Create Convenience Part");

        String jsonString = "{'username':'jsmith','objectName':'Test Convenience Device Part','sequence':'tccctatcagtgatagagattgacatccctatcagtgatagagatactgagcac', 'role':'GENE'}";
        URL url = new URL(this.url + "/create/conveniencePart/");

        String partIDs = HTTPReq(url, jsonString, "POST");
        
        jsonString = "{'username':'jsmith','objectName':'Test Convenience Device','sequence':'actacttcgcatcatgttcatca', 'role':'GENE', 'partIDs':'" + partIDs +"'}";
        url = new URL(this.url + "/create/convenienceDevice/");

        String output = HTTPReq(url, jsonString, "POST");

        System.out.println(output);
    }
//    @Test
//    public void timeToBulkCreate() throws MalformedURLException, IOException, ProtocolException, NoSuchAlgorithmException, KeyManagementException
//    {
//        System.out.println("Testing Bulk Create");
//        URL url = new URL(this.url + "/create/sequence");
//        long start = System.currentTimeMillis();
//        for (int i = 0; i < 5000; i++)
//        {    
//            String jsonString = "{'username':'jsmith','password':'asdf','objectName':'K249" + i + " Sequence','sequence':'atgcagatttatgaaggcaaactgaccgcggaaggcctgcgctttggcattgtggcgagccgctttaaccatgcgc"
//				+ "tggtggatcgcctggtggaaggcgcgattgattgcattgtgcgccatggtggtcgcgaagaagatattaccctggtgcgcgtgccgggcagctgggaaattccggtgg"
//				+ "cggcgggcgaactggcgcgcaaagaagatattgatgcggtgattgcgattggcgtgctgattgaaggcgcggaaccgcattttgattatattgcgagcgaagtgagca"
//				+ "aaggcctggcgaacctgagcctggaactgcgcaaaccgattacctttggcgtgattaccgcggatgaactggaagaagcgattgaacgcgcgggcaccaaacatggca"
//				+ "acaaaggctgggaagcggcgctgagcgcgattgaaatggcgaacctgtttaaaagcctgcgctag'}";
//            
//            HTTPReq(url, jsonString, "POST");
//        }
//        long end = System.currentTimeMillis();
//        System.out.println("Bulk Create in Rest API took " + (end - start) + " MilliSeconds");      
//    }
}
