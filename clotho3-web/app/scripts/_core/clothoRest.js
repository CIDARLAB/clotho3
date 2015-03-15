"use strict";

angular.module('clotho.core')
/**
 * @name ClothoREST
 *
 * @description
 * Service providing REST API access using message formats from ClothoAPI and Socket.
 *
 * todo - expose url in config / globally
 */
  .service('ClothoREST', function ($http, ClothoAuth, Debug) {

    var urlRest = '';
    var Debugger = new Debug('ClothoREST', '#bbbb55');

    var parseForRest = function (msgObject) {
      //todo - parse headers for rest, determine method, format for $http

      return {
        method: 'GET',
        url: urlRest,
        headers: {

        }
      };
    };

    var send = function (msgObject) {
      $http(parseForRest(msgObject));
    };

    return {
      send: send
    };

  });
