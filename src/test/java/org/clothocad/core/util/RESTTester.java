/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import com.mashape.unirest.http.HttpResponse;
import com.mashape.unirest.http.Unirest;
import com.mashape.unirest.http.exceptions.UnirestException;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.InetSocketAddress;
import java.net.MalformedURLException;
import java.net.ProtocolException;
import java.net.Socket;
import java.net.URL;
import java.net.UnknownHostException;
import java.security.KeyManagementException;
import java.security.KeyStoreException;
import java.security.NoSuchAlgorithmException;
import java.security.UnrecoverableKeyException;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.net.ssl.HttpsURLConnection;
import javax.net.ssl.SSLContext;
import javax.net.ssl.TrustManager;
import javax.net.ssl.X509TrustManager;
import org.apache.http.client.HttpClient;
import org.apache.http.conn.ClientConnectionManager;
import org.apache.http.conn.ConnectTimeoutException;
import org.apache.http.conn.scheme.PlainSocketFactory;
import org.apache.http.conn.scheme.Scheme;
import org.apache.http.conn.scheme.SchemeRegistry;
import org.apache.http.conn.scheme.SchemeSocketFactory;
import org.apache.http.impl.client.DefaultHttpClient;
import org.apache.http.impl.conn.SingleClientConnManager;
import org.apache.http.params.HttpParams;
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
        conn.setDoOutput(true);
        conn.setRequestMethod(verb);
        conn.setRequestProperty("Content-Type", "application/json");
        conn.setRequestProperty("Cache-Control", "no-cache");
        conn.setInstanceFollowRedirects(false);

        OutputStream os = conn.getOutputStream();
        os.write(jsonString.getBytes());
        os.flush();

        if (conn.getResponseCode() == 200) {
            System.out.println("SUCCESS!");

            //print result
            BufferedReader br = new BufferedReader(new InputStreamReader((conn.getInputStream())));

            String output;
            String alloutput = "";
            while ((output = br.readLine()) != null) {
                alloutput += output;
            }
            conn.disconnect();
            return alloutput;
        }
        conn.disconnect();
        return "ERROR";
    }
            

    @Test
    public void testCreateUser() throws MalformedURLException, IOException, KeyManagementException, NoSuchAlgorithmException {
        System.out.println("Testing Create User");    
        String jsonString = "{'username':'jsmith','password':'asdf'}";
        URL url = new URL(this.url + "/create/user");
       
        String output = HTTPReq(url, jsonString, "POST");
        
        System.out.println(output);
    }
    
    @Test
    public void testCreateSequence() throws MalformedURLException, IOException, KeyManagementException, NoSuchAlgorithmException {
        System.out.println("Testing Create Sequence");
        String jsonString = "{'username':'jsmith','password':'asdf','objectName':'Test Sequence','sequence':'ata'}";
        URL url = new URL(this.url + "/create/sequence");
        
        String output = HTTPReq(url, jsonString, "POST");
        
        System.out.println(output);
    }
    
    @Test
<<<<<<< HEAD
    public void testCreatePart() throws MalformedURLException, IOException, KeyManagementException, NoSuchAlgorithmException {
        
        String jsonString = "{'username':'jsmith','password':'asdf','objectName':'Test Sequence','sequence':'ata'}";
        URL url = new URL(this.url + "/create/sequence");      
        String seqId = HTTPReq(url, jsonString, "POST");
               
        System.out.println("Testing Create Part");
        jsonString = "{'username':'jsmith','password':'asdf','objectName':'Test Part', 'id':'"+ seqId +"'}";
        url = new URL(this.url + "/create/part");
        
        String output = HTTPReq(url, jsonString, "POST");
        
=======
    public void testGetByName() throws MalformedURLException, IOException, KeyManagementException, NoSuchAlgorithmException {
        System.out.println("Testing Get Sequence by Name");
        
        String jsonString = "{'username':'jsmith','password':'asdf','objectName':'Test Sequence','sequence':'ata'}";
        URL url = new URL(this.url + "/create/sequence");
        
        
        jsonString = "{'username':'jsmith','password':'asdf','objectName':'Test Sequence'}";
        url = new URL(this.url + "/get/getByName");
        
        String output = HTTPReq(url, jsonString, "GET");
                
>>>>>>> 7e9ae75ef19e53448011d7384524d83111ad9510
        System.out.println(output);
    }
    
    @Test
<<<<<<< HEAD
    public void testGetById() throws MalformedURLException, IOException, KeyManagementException, NoSuchAlgorithmException {      
        
        String jsonString = "{'username':'jsmith','password':'asdf','objectName':'Test Sequence','sequence':'ata'}";
        URL url = new URL(this.url + "/create/sequence");      
        String seqId = HTTPReq(url, jsonString, "POST");
        
        System.out.println("Testing Get By Id");
        jsonString = "{'username':'jsmith','password':'asdf','id':'"+ seqId +"'}";
        url = new URL(this.url + "/get/getById");
        
        String output = HTTPReq(url, jsonString, "GET");
=======
    public void testDelete() throws MalformedURLException, IOException, KeyManagementException, NoSuchAlgorithmException {
        System.out.println("Testing Delete Sequence");
        
        String jsonString = "{'username':'jsmith','password':'asdf','objectName':'Test Sequence','sequence':'ata'}";
        URL url = new URL(this.url + "/create/sequence");
        
        String sequenceId = HTTPReq(url, jsonString, "POST");;
        
        
        jsonString = "{'username':'jsmith','password':'asdf','id':" + sequenceId + "}";
        url = new URL(this.url + "/delete/delete");
        
        String output = HTTPReq(url, jsonString, "DELETE");
>>>>>>> 7e9ae75ef19e53448011d7384524d83111ad9510
        
        System.out.println(output);
    }
}
