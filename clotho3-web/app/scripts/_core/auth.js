
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
			getCurrentUserDeferred = [];

		//listen for login and logout events, update info and resolve requests
		PubSub.on('auth:login auth:logout', function handleLoginEvent (info) {
			//null or undefined suggest logout()
			Debugger.log('login event info:', info);
			if (angular.isEmpty(info)) {
				currentUser = null;
			} else {
				currentUser = info;
				//process the deferred queue only on login()
				while (getCurrentUserDeferred.length > 0) {
					var def = getCurrentUserDeferred.pop();
					def.resolve(currentUser);
				}
			}
		});

		PubSub.on('auth:error', function handleLoginError () {

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
			}
		};
	});