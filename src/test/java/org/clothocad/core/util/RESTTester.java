/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import com.mashape.unirest.http.HttpResponse;
import com.mashape.unirest.http.Unirest;
import com.mashape.unirest.http.exceptions.UnirestException;
import java.io.IOException;
import java.net.InetSocketAddress;
import java.net.Socket;
import java.net.UnknownHostException;
import java.security.KeyManagementException;
import java.security.KeyStoreException;
import java.security.NoSuchAlgorithmException;
import java.security.UnrecoverableKeyException;
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
import org.junit.Test;

/**
 *
 * @author David
 */
public class RESTTester {
    
    @Test
    public void testCreateUser() throws UnirestException{
        String url = "https://localhost:8080/data/post/createUser";
//      "{"username":"jasmith","password":"asdf","type":"sequence","name":"B34 Sequence","value":"atag"}";
    
        HttpResponse<String> res = Unirest.post(url)
                .header("accept", "application/json")
                .field("username", "jasmith")
                .field("password", "asdf")
                .field("type","sequence")
                .field("name", "B34 Sequence")
                .field("value", "atag").asString();
        
        String body = res.getBody();
        
        
        System.out.println(body);
    }
    
}
