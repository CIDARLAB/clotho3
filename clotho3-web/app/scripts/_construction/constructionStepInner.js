angular.module('clotho.construction')
/**
 * @ngdoc directive
 * @name constructionStepInner
 *
 * @restrict A
 *
 * @description
 * Grabs template specific to the construction file step. Does not create isolate scope. Pass reaction name as attr.
 *
 * @example
 * <div construction-step-inner="{{step.reaction}}"></div>
 */
	.directive('constructionStepInner', function ($http, $compile, ConstructionReactions) {
		return {
			restrict : 'A',
			link : function (scope, element, attrs) {
				scope.$watch(function () {
					return attrs.constructionStepInner
				}, function (val) {
					if (!!val) {
						var reaction = ConstructionReactions[val];

						!angular.isEmpty(reaction) && $http.get(reaction.template, {cache : true})
						.then(function (data) {
							element.html($compile(data.data)(scope));
						});
					}
				});
			}
		}
	});