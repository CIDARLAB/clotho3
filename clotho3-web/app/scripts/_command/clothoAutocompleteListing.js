angular.module('clotho.tokenizer')
	/*
	 * internal directive which displays the actual list of autocompletions
	 *
	 * replies on each object in list having a UNQIUE 'id' field
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
					if (scope.forceVisible === true || scope.forceVisible === false) {
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