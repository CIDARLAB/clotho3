
angular.module('clotho.core')
/**
 * @name ClothoAuth
 *
 * @description
 * Singleton to handle user authentication and credentials
 *
 * todo - secure credentials and user info
 * todo - expose on $rootScope - see firebaseSimpleLogin
 *
 */
	.service('ClothoAuth', function(PubSub, Debug, $q) {

		var Debugger = new Debug('ClothoAuth', '#EE3333');

		var userId = null,                      //string
      credentials = null,                   //object
      stateListeners = [],
			getCurrentUserDeferred = [];

    function cloneCurrentUser () {
      return userId;
    }

    function cloneCredentials() {
      return angular.copy(credentials);
    }

    function triggerStateListeners (state) {
      angular.forEach(stateListeners, function (callback) {
        callback.call(null, state, cloneCurrentUser());
      });
    }

		PubSub.on('auth:login', function handleLoginEvent (info) {
			Debugger.log('login event info:', info);

      //todo - update once server responds correctly
      userId = info.primaryEmail;
      credentials = info.credentials;

      //process the deferred queue only on login()
      while (getCurrentUserDeferred.length > 0) {
        var def = getCurrentUserDeferred.pop();
        def.resolve(cloneCurrentUser());
      }

      triggerStateListeners('login');
		});

    PubSub.on('auth:logout', function handleLogoutEvent (info) {
      Debugger.log('logout', info);
      userId = null;
      credentials = null;
      triggerStateListeners('logout');
    });

		PubSub.on('auth:error', function handleLoginError () {
      triggerStateListeners('error');
		});

		return {
      isLoggedIn: function() {
				return !angular.isEmpty(userId);
			},

			getCurrentUser: cloneCurrentUser,

      getCredentials: cloneCredentials,

      /**
       * @name onLogin
       *
       * @description
       * Returns a promise, which is fulfilled on login() with userId
       *
       * If not logged in, will register promise and wait until login event. This means the promise may never be fulfilled.
       *
       * @returns {Promise} Fulfilled with userId
       */
      onLogin: function () {
        if (angular.isDefined(userId)) {
          $q.when(cloneCurrentUser());
        } else {
          var deferred = $q.defer();
          getCurrentUserDeferred.push(deferred);
          return deferred.promise;
        }
      },

      /**
       * @name addStateListener
       * @description
       * Register a callback when state changes. To listen to a particular state, use PubSub and listen for events 'auth:login', 'auth:logout', 'auth:error'.
       *
       * @param callback {Function}
       * passed state {string},
       *
       * @returns {Function} deregistration function
       */
      addStateListener : function (callback) {
        //PubSub handles check for callback is function
        return PubSub.on('auth:login auth:logout auth:error', callback);
      }
		};
	});
