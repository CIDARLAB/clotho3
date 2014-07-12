package org.clothocad.core.testers;

import javax.jms.*;
import org.fusesource.stomp.jms.StompJmsConnectionFactory;

public class Producer {

    private static final String BROKER_URL = "tcp://localhost:61613";
    private static final Boolean NON_TRANSACTED = false;
    private static final int NUM_MESSAGES_TO_SEND = 100;
    private static final long DELAY = 100;

    public static void main(String[] args) {

        StompJmsConnectionFactory factory =
                new StompJmsConnectionFactory();
        factory.setBrokerURI(BROKER_URL);
        Connection connection = null;

        try {

            connection = factory.createConnection("admin", "password");
            connection.start();

            System.out.println("Started Connection...");
            Session session = connection.createSession(NON_TRANSACTED, Session.AUTO_ACKNOWLEDGE);
            Destination destination = session.createQueue("CLOTHO");
            MessageProducer producer = session.createProducer(destination);

            for (int i = 0; i < NUM_MESSAGES_TO_SEND; i++) {
                TextMessage message = session.createTextMessage("Message #" + i);
                System.out.println("Sending message #" + i);
                producer.send(message);
                Thread.sleep(DELAY);
            }

            producer.close();
            session.close();

        } catch (Exception e) {
            System.out.println("Caught exception!");
        } finally {
            if (connection != null) {
                try {
                    connection.close();
                } catch (JMSException e) {
                    System.out.println("Could not close an open connection...");
                }
            }
        }
    }
}
