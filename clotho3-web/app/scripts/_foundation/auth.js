'use strict';

angular.module('clotho.foundation')
/**
 * A directive that shows elements only when the given authentication state is in effect
 *
 * <code>
 *    <div clotho-show-auth="login">{{auth.user.id}} is logged in</div>
 *    <div clotho-show-auth="logout">Logged out</div>
 *    <div clotho-show-auth="error">An error occurred: {{auth.error}}</div>
 *    <div clotho-show-auth="logout,error">This appears for logout or for error condition!</div>
 * </code>
 */
.directive('clothoShowAuth', function (PubSub) {
	var loginState = 'logout';
	PubSub.on('auth:login',  function() { loginState = 'login'; });
	PubSub.on('auth:logout', function() { loginState = 'logout'; });
	PubSub.on('auth:error',  function() { loginState = 'error'; });

	function getExpectedState(scope, attr) {
		var expState = scope.$eval(attr);
		if( typeof(expState) !== 'string' && !angular.isArray(expState) ) {
			expState = attr;
		}
		if( typeof(expState) === 'string' ) {
			expState = expState.split(',');
		}
		return expState;
	}

	function inList(needle, list) {
		var res = false;
		angular.forEach(list, function(x) {
			if( x === needle ) {
				res = true;
				return true;
			}
			return false;
		});
		return res;
	}

	function assertValidStates(states) {
		if( !states.length ) {
			throw new Error('clotho-show-auth directive must be login, logout, or error (you may use a comma-separated list)');
		}
		angular.forEach(states, function(s) {
			if( !inList(s, ['login', 'logout', 'error']) ) {
				throw new Error('Invalid state "'+s+'" for clotho-show-auth directive, must be one of login, logout, or error');
			}
		});
		return true;
	}

	return {
		restrict: 'A',
		link: function(scope, el, attr) {
			var expState = getExpectedState(scope, attr.clothoShowAuth);
			assertValidStates(expState);
			function fn() {
				var show = inList(loginState, expState);
				// sometimes if ngCloak exists on same element, they argue, so make sure that
				// this one always runs last for reliability
				setTimeout(function() {
					el.toggleClass('ng-cloak', !show);
				}, 0);
			}
			fn();
			PubSub.on('auth:login auth:logout auth:error',  fn);
		}
	};
});