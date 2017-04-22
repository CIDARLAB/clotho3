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
import org.json.JSONObject;
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

//    //@Test
//    public void testCreateUser() throws MalformedURLException, IOException, KeyManagementException, NoSuchAlgorithmException {
//        System.out.println("Testing Create User");
//        String jsonString = "{'username':'jsmith400','password':'asdf'}";
//        URL url = new URL(this.url + "/create/user");
//
//        String output = HTTPReq(url, jsonString, "POST");
//
//        System.out.println(output);
//    }
//    @Test
    public void testCreateSequence() throws MalformedURLException, IOException, KeyManagementException, NoSuchAlgorithmException {
        System.out.println("Testing Create Sequence");
        String jsonString = "{'username':'jsmith', 'objectName':'Test Sequence','description':'Test Sequence', 'sequence':'ata'}";
        URL url = new URL(this.url + "/create/sequence/");

        String output = HTTPReq(url, jsonString, "POST");

        System.out.println(output);
    }

    //@Test
    public void testCreatePart() throws MalformedURLException, IOException, KeyManagementException, NoSuchAlgorithmException {

        String jsonString = "{'username':'jsmith','objectName':'Test Sequence','sequence':'ata'}";
        URL url = new URL(this.url + "/create/sequence");
        String seqId = HTTPReq(url, jsonString, "POST");

        System.out.println("Testing Create Part");
        jsonString = "{'username':'jsmith','objectName':'Test Part', 'description':'Test Part', 'id':'" + seqId + "'}";
        url = new URL(this.url + "/create/part");

        String output = HTTPReq(url, jsonString, "POST");

        System.out.println(output);
    }

//    @Test
    public void testGetByName() throws MalformedURLException, IOException, KeyManagementException, NoSuchAlgorithmException {
        System.out.println("Testing Get Sequence by Name");

        String jsonString = "{'username':'jsmith','objectName':'Test Sequence', 'description':'Test Sequence', 'sequence':'ata'}";
        URL url = new URL(this.url + "/create/sequence");
        String seqId = HTTPReq(url, jsonString, "POST");

        url = new URL("https://localhost:8443/data/get/getByName/Test%20Sequence///20");
        String output = HTTPReq(url, "", "GET");
        JSONObject obj = new JSONObject(output);
        System.out.println(output);
        Object pages = obj.get("page_count");
        while (obj.getJSONArray("links").getJSONObject(0).has("next")) {
            url = new URL("https://localhost:8443/data/get/getByName/Test%20Sequence" + obj.getJSONArray("links").getJSONObject(0).getString("next"));
            output = HTTPReq(url, "", "GET");
            obj = new JSONObject(output);
            System.out.println(output);
        }
    }

//    @Test
    public void testGetById() throws MalformedURLException, IOException, KeyManagementException, NoSuchAlgorithmException {
        System.out.println("Testing Get By Id");

        String jsonString = "{'username':'jsmith','objectName':'Test Sequence','sequence':'ata'}";
        URL url = new URL(this.url + "/create/sequence");
        String seqId = HTTPReq(url, jsonString, "POST");

        url = new URL("https://localhost:8443/data/get/getByID/" + seqId);

        String output = HTTPReq(url, "", "GET");

        System.out.println(output);

    }

    //@Test
    public void testSet() throws MalformedURLException, IOException, KeyManagementException, NoSuchAlgorithmException {

        String jsonString = "{'username':'jsmith','objectName':'Test Sequence','sequence':'ata'}";
        URL url = new URL(this.url + "/create/sequence");
        String seqId = HTTPReq(url, jsonString, "POST");

        System.out.println("Testing Set");

        jsonString = "{'username':'jsmith','sequence' : 'atatatatatatat','id' : '" + seqId + "'}";
        url = new URL(this.url + "/set");

        String output = HTTPReq(url, jsonString, "PUT");

        System.out.println(output);
    }

//    @Test
    public void testDelete() throws MalformedURLException, IOException, KeyManagementException, NoSuchAlgorithmException {
        System.out.println("Testing Delete Sequence");

        String jsonString = "{'username':'jsmith','objectName':'Test Sequence','sequence':'ata'}";
        URL url = new URL(this.url + "/create/sequence");

        String sequenceId = HTTPReq(url, jsonString, "POST");

        url = new URL(this.url + "/delete/" + sequenceId);

        String output = HTTPReq(url, jsonString, "DELETE");

        System.out.println(output);
    }

//    @Test
    public void testConveniencePart() throws MalformedURLException, IOException, KeyManagementException, NoSuchAlgorithmException {
        System.out.println("Testing Create Convenience Part");

//        String jsonString = "{'username':'jsmith','objectName':'Test Convenience Part','displayID':'Test Convenience Part','sequence':'tccctatcagtgatagagattgacatccctatcagtgatagagatactgagcac', 'role':'GENE', 'params': [{'name':'n', 'value':'121.5', 'variable':'var', 'units' : 'unit'}]}";
        String jsonString = "{'username':'jsmith','objectName':'Test Convenience Part'}";

        URL url = new URL(this.url + "/create/conveniencePart/");

        String output = HTTPReq(url, jsonString, "POST");

        System.out.println(output);

//        System.out.println("Testing Get Convenience Part");
//
//        url = new URL(this.url + "/get/conveniencePart/");
//        output = HTTPReq(url, jsonString, "POST");
//
//        System.out.println(output);
    }

    @Test
    public void testConvenienceDevice() throws MalformedURLException, IOException, KeyManagementException, NoSuchAlgorithmException {
        System.out.println("Testing Create Convenience Device");

        String jsonString1 = "{'username':'jsmith','objectName':'Test Convenience Device Part1','displayID':'Test Convenience Device Part1','sequence':'tccctatcagtgatagagattgacatccctatcagtgatagagatactgagcac', 'role':'GENE', 'params': [{'name':'n', 'value':'121.5', 'variable':'var', 'units' : 'unit'}]}";
        URL url1 = new URL(this.url + "/create/conveniencePart/");

        String jsonString2 = "{'username':'jsmith','objectName':'Test Convenience Device Part2','displayID':'Test Convenience Device Part2', 'sequence':'tccctatcagtgatagagattgacatccctatcgagatactgagcac', 'role':'GENE', 'params': [{'name':'n', 'value':'121.5', 'variable':'var', 'units' : 'unit'}]}";
        URL url2 = new URL(this.url + "/create/conveniencePart/");

        String partID1 = HTTPReq(url1, jsonString1, "POST");
        String partID2 = HTTPReq(url2, jsonString2, "POST");

        String partIDs = partID1 + "," + partID2;

        String jsonString = "{'username':'jsmith','objectName':'Test Convenience Device','displayID':'Test Convenience Device', 'createSeqFromParts':'False', 'partIDs':'" + partIDs + "'}";
        URL url = new URL(this.url + "/create/convenienceDevice/");

        String output = HTTPReq(url, jsonString, "POST");

        System.out.println(output);

//        System.out.println("Testing Get Convenience Device");
//
//        jsonString = "{'username':'jsmith','objectName':'Test Convenience Device','displayID':'Test Convenience Device', 'sequence':'tccctatcagtgatagagattgacatccctatcagtgatagagatactgagcac','role':'GENE','params':[{'name':'n', 'value':'121.5', 'variable':'var', 'units' : 'unit'}], 'parts':[{'name':'Test Convenience Device Part1', 'description':'Test Convenience Device Part1'}]}";
//        url = new URL(this.url + "/get/convenienceDevice/");
//        output = HTTPReq(url, jsonString, "POST");
//
//        System.out.println(output);
    }

//    //@Test
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