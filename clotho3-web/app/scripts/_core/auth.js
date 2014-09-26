
angular.module('clotho.core')
/**
 * @name ClothoAuth
 *
 * @description
 * Singleton to handle user authentication and credentials
 *
 * todo - this is not secure, because anyone could change the credentials. Shouldn't be abe to do anything on the server, but still...
 *
 */
	.service('ClothoAuth', function(PubSub, Debug, $q) {

		var Debugger = new Debug('ClothoAuth', '#EE3333');

		var currentUser = null,
      stateListeners = [],
			getCurrentUserDeferred = [];

    function triggerStateListeners (state) {
      angular.forEach(stateListeners, function (callback) {
        callback.call(null, state, currentUser);
      });
    }

		PubSub.on('auth:login', function handleLoginEvent (info) {
			Debugger.log('login event info:', info);
      currentUser = info;
      //process the deferred queue only on login()
      while (getCurrentUserDeferred.length > 0) {
        var def = getCurrentUserDeferred.pop();
        def.resolve(currentUser);
      }
      triggerStateListeners('login');
		});

    PubSub.on('auth:logout', function handleLogoutEvent (info) {
      Debugger.log('logout', info);
      currentUser = null;
      triggerStateListeners('logout');
    });

		PubSub.on('auth:error', function handleLoginError () {
      triggerStateListeners('error');
		});

		return {
      //for data binding
      //NOT SECURE AT ALL
      currentUser : currentUser,

			isLoggedIn: function() {
				return !angular.isEmpty(currentUser);
			},

			// Gets a promise for the current user info.
			getCurrentUser: function() {
				if (angular.isDefined(currentUser)) {
					$q.when(currentUser);
				} else {
					var deferred = $q.defer();
					getCurrentUserDeferred.push(deferred);
					return deferred.promise;
				}
			},

      addStateListener : function (callback) {
        angular.isFunction(callback) && stateListeners.push(callback)
      }
		};
	});
