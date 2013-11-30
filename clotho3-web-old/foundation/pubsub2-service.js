'use strict';

Application.Foundation.service('PubSub', function() {
    /*TODO - rewrite with following goals:

     - Track Multiple callbacks for each event
     - Track reference for each callback
     - On() Pass multiple events, space-separated
     - Once() via Flag for Boolean 'once'
     - Each() (unexposed) to handle multiple/single items the same way
     - Appropriate handle return
     - Off() by handle (event, specific callback)
     - Destroy() by reference
     - Clear() by event

     Note: What this means
     - two way binding between: event, reference, callback

     */

    //see if already exists, don't re-instantiate
    return (window.$clotho.$pubsub2) ? window.$clotho.$pubsub2 : window.$clotho.$pubsub2 = generatePubSubObject();

    function generatePubSubObject() {
        /*****
         CONFIG
         *******/

        var map = {};

        //split events passed in by a space
        var eventSplitter = /\s+/;

        /*********
         Functions
         **********/

        var logListeners = function() {
            console.log(map);
        };

        var trigger = function(topic, args) {

        };

        var on = function(topic, callback, ref, one) {

        };

        var once = function(topic, callback, ref) {
            on(topic, callback, ref);
        };

        var off = function() {

        };

        var destroy = function(ref) {

        };

        var clear = function(topic) {

        };

        return {
            logListeners : logListeners,
            trigger: trigger,
            on : on,
            once : once,
            off : off,
            destroy : destroy,
            clear : clear
        }
    }
});