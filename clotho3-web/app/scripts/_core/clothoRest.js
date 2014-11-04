"use strict";

angular.module('clotho.core')
/**
 * @name ClothoREST
 *
 * @description
 * Service providing REST API access using message formats from ClothoAPI and Socket.
 *
 * See https://docs.angularjs.org/api/ng/service/$http for usage
 *
 */
.provider('ClothoREST', function ClothoRESTProvider() {
  var self = this;
  self.rootUrl = 'https://localhost:8443/data/';
  self.defaults = {
    cache: false,
    timeout: 5000
  };


  this.$get = function ClothoRESTFactory($http, ClothoAuth, Debug, Base64Util) {

    var Debugger = new Debug('ClothoREST', '#bbbb55');

    var prefixUrl = function (endpoint) {
      return self.rootUrl + (angular.isString(endpoint) ? endpoint : '');
    };

    var createHeaders = function (config) {
      var userId = ClothoAuth.getUserId();
      if (!userId) {
        //not logged in
        throw "Not authorized - you are not logged in";
      }

      //todo - update format
      var headers = {
        'Authentication' : Base64Util.encode(userId + ":" + ClothoAuth.getToken())
      };

      var params = {
        headers: headers
      };

      return angular.extend({}, self.defaults, params, config);
    };

    /**
     * Directly access $http to make a request. You must specify URL and method yourself
     * @param config
     * @returns {Promise}
     */
    var send = function (config) {
      return $http(createHeaders(config));
    };

    var rest_create = function (data, config) {
      return send(angular.extend({
        method: "POST",
        url: prefixUrl(),
        data: data
      }, config))
    };

    var rest_get = function (id, config) {
      return send(angular.extend({
        method: 'GET',
        url : prefixUrl(id)
      }, config));
    };

    /**
     *
     * @param queryParams {Object|String}
     * @param config {Object}
     * @returns {Promise}
     */
    var rest_query = function (queryParams, config) {
      return send(angular.extend({
        method: 'GET',
        url : prefixUrl(),
        params: queryParams
      }, config));
    };

    var rest_destroy = function (id, config) {
      return send(angular.extend({
        method: 'DELETE',
        url : prefixUrl(id)
      }, config));
    };

    var rest_set = function (id, data, config) {
      return send(angular.extend({
        method: 'PUT',
        url : prefixUrl(id),
        data : data
      }, config));
    };

    var rest_run = function (id, data, config) {
      return send(angular.extend({
        method: 'PUT',
        url : prefixUrl(id),
        data: data
      }, config));
    };

    return {
      'send' : send,
      'create' : rest_create,
      'get' : rest_get,
      'query': rest_query,
      'destroy' : rest_destroy,
      'set' : rest_set,
      'run' : rest_run
    };
  };
});
