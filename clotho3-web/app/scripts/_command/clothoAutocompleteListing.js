angular.module('clotho.tokenizer')
	/*
	 * internal directive which displays the actual list of autocompletions
	 */
	.directive('clothoAutocompleteListing', function () {
		return {
			restrict:'EA',
			scope:{
				autocompletions:'=',
				query:'=',
				active:'=',
				hasFocus: '=',
				select:'&',
				passedPlacement : '@?',
				forceVisible : '@?'
			},
			replace:true,
			templateUrl:'views/_command/autocompleteListing.html',
			link:function (scope, element, attrs) {

				scope.isVisible = function () {
					if (angular.isBoolean(scope.forceVisible)) {
						return scope.forceVisible;
					} else {
						return scope.hasFocus && scope.autocompletions.length;
					}
				};

				scope.isOpen = function () {
					return scope.hasFocus && scope.matches.length > 0;
				};

				scope.isActive = function (matchIdx) {
					return scope.active == matchIdx;
				};

				scope.selectActive = function (matchIdx) {
					scope.active = matchIdx;
				};

				scope.selectMatch = function (activeIdx) {
					scope.select({activeIdx:activeIdx});
				};
			}
		};
	});