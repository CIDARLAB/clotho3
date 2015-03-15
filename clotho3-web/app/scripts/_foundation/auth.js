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
 *
 * Note - you can use ng-disabled and a scope variable tied to the userInfo to handle dynamically disabling form elements.. doesn't really make sense to create a directive to handle that directly
 */
.directive('clothoShowAuth', function (PubSub, ClothoAuth) {

  var validstates = ['login', 'logout', 'error'];

  //todo - may not want to eval if just passing string, not scope variables. add $observe
	function getExpectedStates(scope, attr) {
		var expState = scope.$eval(attr);
		if( !angular.isString(expState) && angular.isArray(expState) ) {
			expState = attr;
		}
		if( angular.isString(expState) ) {
			expState = expState.split(',');
		}
		return expState;
	}

	function inList(needle, list) {
    for (var i = 0; i < list.length; i++) {
      if (list[i] === needle) {
        return true;
      }
    }
    return false;
	}

	function assertValidStates(states) {
		if( !states.length ) {
			throw new Error('clotho-show-auth directive must be (you may use a comma-separated list): ' + validstates.join(', '));
		}
		angular.forEach(states, function(s) {
			if( !inList(s, validstates) ) {
				throw new Error('Invalid state(s) "'+s+'" for clotho-show-auth directive, must be one of: ' + validstates.join(', '));
			}
		});
		return true;
	}

	return {
		restrict: 'A',
		link: function(scope, el, attr) {
			var expStates = getExpectedStates(scope, attr.clothoShowAuth);
			assertValidStates(expStates);
			function fn(newState) {
				var show = inList(newState, expStates);
				// sometimes if ngCloak exists on same element, they argue, so make sure that
				// this one always runs last for reliability
				setTimeout(function() {
					el.toggleClass('ng-cloak', !show);
				}, 0);
			}
      //start in logout state
			//fn('logout');
      ClothoAuth.addStateListener(fn);
		}
	};
});
