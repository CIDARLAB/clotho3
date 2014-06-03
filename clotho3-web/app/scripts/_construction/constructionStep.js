angular.module('clotho.construction')
/**
 * @ngdoc directive
 * @name constructionStep
 *
 * @restrict A
 *
 * @description
 * Visualize a step of a construction file. Used internally in constructionFileView directive.
 *
 * @note
 * expects 'step' on scope, and $index as attr value
 *
 * @example
 * <div ng-repeat="step in file.steps"
        construction-step></div>
 */
	.directive('constructionStep', function (ConstructionReactions, $parse) {
		return {
			restrict : 'A',
			templateUrl : 'views/_construction/constructionStep.html',
			link : function (scope, element, attrs) {
				scope.$watch('step', function (newval) {
					scope.reaction = ConstructionReactions[newval.reaction];
				});
			}
		}
	})
/**
 * @name constructionStepInner
 * *
 * @description
 * Used internally by constructionStep directive
 * Grabs template specific to the construction file step. Does not create isolate scope. Pass reaction name as attr.
 *
 * @example
 * (given 'reaction' on parent scope)
 * <div construction-step-inner></div>
 */
	.directive('constructionStepInner', function ($http, $compile, ConstructionReactions) {
		return {
			restrict : 'A',
			link : function (scope, element, attrs) {
				scope.$watch('reaction.template', function (newval) {
					angular.isDefined(newval) && $http.get(newval, {cache : true})
					.then(function (data) {
						element.html($compile(data.data)(scope));
					});
				});
			}
		}
	});