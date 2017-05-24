package org.clothocad.core;

import org.clothocad.core.persistence.DBClassLoader;
import org.clothocad.core.persistence.IdUtils;
import org.clothocad.core.schema.Schema;

import com.google.inject.AbstractModule;
import com.google.inject.name.Names;
import com.google.inject.Provides;
import java.math.BigInteger;
import java.security.KeyPair;
import java.security.KeyPairGenerator;
import java.security.KeyStore;
import java.security.SecureRandom;
import java.security.Security;
import java.security.cert.Certificate;
import java.security.cert.X509Certificate;
import java.util.Date;


import org.eclipse.jetty.util.ssl.SslContextFactory;

import java.util.Properties;

import javax.inject.Singleton;
import javax.security.auth.x500.X500Principal;
import org.bouncycastle.jce.provider.BouncyCastleProvider;
import org.bouncycastle.x509.X509V1CertificateGenerator;

public class ClothoModule extends AbstractModule {

    protected final Properties config;

    public ClothoModule(Properties config) {
        this.config = config == null ?
            ConfigOption.getDefaultConfig() : config;
    }

    @Override
    protected void configure() {
        Names.bindProperties(binder(), config);
        requestStaticInjection(Schema.class);
        requestStaticInjection(IdUtils.class);
        bind(DBClassLoader.class).in(Singleton.class);
    }

    @Provides
    protected SslContextFactory provideSslConnector() throws Exception {
        SslContextFactory cf = new SslContextFactory();
        
//        cf.setKeyStorePath(config.getProperty("keystorepath"));
//        cf.setKeyStorePassword(config.getProperty("keystorepass"));
        
        cf.setKeyStore(provideKeyStore());



        //SslSelectChannelConnector sslConnector = new SslSelectChannelConnector(cf);
        return cf;
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
