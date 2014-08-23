
angular.module('clotho.core')
/**
 * @name ClothoAuth
 *
 * @description
 * Singleton to handle user authentication and credentials
 *
 */
	.service('ClothoAuth', function(Clotho, PubSub, $q) {

		//todo - store credentials once server returns them - GH#394

		var currentUser = null,
			waitForAuth = $q.defer(),
			getCurrentUserDeferred = [];

		//listen for login and logout events, update info and resolve requests
		PubSub.on('auth:login auth:logout', handleLoginEvent);
		function handleLoginEvent (info) {
			//take the info... if it's null or undefined then logout
			//todo - check for error
			currentUser = angular.isDefined(info) ? angular.copy(info) : null;

			//process the deferred queue
			while (getCurrentUserDeferred.length > 0) {
				var def = getCurrentUserDeferred.pop();
				def.resolve(currentUser);
			}
		}

		return {

			login: Clotho.login,
			logout: Clotho.logout,

			isLoggedIn: function() {
				return angular.isDefined(currentUser);
			},

			// Gets a promise for the current user info.
			getCurrentUser: function() {
				var deferred = $q.defer();

				if (angular.isDefined(currentUser)) {
					deferred.resolve(currentUser);
				} else {
					getCurrentUserDeferred.push(deferred);
				}

				return deferred.promise;
			}
		};
	});