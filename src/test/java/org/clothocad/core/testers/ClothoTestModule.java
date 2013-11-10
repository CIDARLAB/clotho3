/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.testers;

import java.math.BigInteger;
import java.security.KeyPair;
import java.security.KeyPairGenerator;
import java.security.KeyStore;
import java.security.SecureRandom;
import java.security.Security;
import java.security.cert.Certificate;
import java.security.cert.X509Certificate;
import java.util.Date;
import java.util.Properties;
import javax.security.auth.x500.X500Principal;
import org.bouncycastle.jce.provider.BouncyCastleProvider;
import org.bouncycastle.x509.X509V1CertificateGenerator;
import org.clothocad.core.ClothoModule;
import org.clothocad.core.communication.Router;
import org.clothocad.core.util.TestRouter;
import org.eclipse.jetty.server.ssl.SslConnector;
import org.eclipse.jetty.server.ssl.SslSelectChannelConnector;
import org.eclipse.jetty.util.ssl.SslContextFactory;

/**
 *
 * @author spaige
 */
public class ClothoTestModule extends ClothoModule {

    public ClothoTestModule(Properties props) {
        super(props);
        defaults.put("dbname", "testClotho");
    }

    public ClothoTestModule() {
        this(null);
    }

    @Override
    protected void configure() {
        super.configure(); //To change body of generated methods, choose Tools | Templates.
        bind(Router.class).to(TestRouter.class);
    }

    @Override
    protected SslConnector provideSslConnector() throws Exception {
        SslContextFactory cf = new SslContextFactory();
        cf.setKeyStore(provideKeyStore());
        SslSelectChannelConnector sslConnector = new SslSelectChannelConnector(cf);                
        return sslConnector;
    }
    
    protected KeyStore provideKeyStore() throws Exception {
        Security.addProvider(new BouncyCastleProvider());

        //provides keystore with self-signed certificate that is never saved to disk
        //source: http://www.mayrhofer.eu.org/create-x509-certs-in-java
        SecureRandom random = new SecureRandom();
        //Jetty uses the SunJSSE algorithms
        String providerName = "SunJSSE";

        //to get key info out of the generated cert, we also have to generate
        //the keypair with BC
        KeyPairGenerator keyGen = KeyPairGenerator.getInstance("RSA", "BC");
        keyGen.initialize(1024, random);
        KeyPair keypair = keyGen.generateKeyPair();


        // GENERATE THE X509 CERTIFICATE (source: http://code.google.com/p/xebia-france/wiki/HowToGenerateaSelfSignedX509CertificateInJava)
        // yesterday
        Date validityBeginDate = new Date(System.currentTimeMillis() - 24 * 60 * 60 * 1000);
        // in 2 years
        Date validityEndDate = new Date(System.currentTimeMillis() + new Long(2) * 365 * 24 * 60 * 60 * 1000);


        X509V1CertificateGenerator certGen = new X509V1CertificateGenerator();
        X500Principal dnName = new X500Principal("CN=localhost");

        certGen.setSerialNumber(BigInteger.valueOf(System.currentTimeMillis())); //XXX?
        certGen.setSubjectDN(dnName);
        certGen.setIssuerDN(dnName);
        certGen.setNotBefore(validityBeginDate);
        certGen.setNotAfter(validityEndDate);
        certGen.setPublicKey(keypair.getPublic());
        certGen.setSignatureAlgorithm("SHA1withRSA");

        X509Certificate cert = certGen.generate(keypair.getPrivate(), providerName);

        cert.checkValidity();
        cert.verify(keypair.getPublic());
        cert.getPublicKey();

        KeyStore keystore = KeyStore.getInstance("PKCS12", providerName);
        keystore.load(null, null);
        //http://www.docjar.com/html/api/org/bouncycastle/jce/examples/PKCS12Example.java.html
        keystore.setKeyEntry("clotho", keypair.getPrivate(), null, new Certificate[]{cert});

        return keystore;

    }
}
