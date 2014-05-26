angular.module('clotho.construction')
/**
 * @ngdoc directive
 * @name constructionStep
 *
 * @restrict A
 *
 * @description
 * Visualize a step of a construction file. Used internally in constructionFileView directive.
 */
	.directive('constructionStep', function (ConstructionReactions, $parse) {
		return {
			restrict : 'A',
			templateUrl : 'views/_construction/constructionStep.html',
			link: function (scope, element, attrs, constructionFileCtrl) {
				scope.$watch(function () {
					return $parse(attrs.constructionStep)(scope)
				}, function (val) {
					if (angular.isNumber(val)) {
						scope.step = scope.file.steps[val];
						scope.reaction = ConstructionReactions[scope.step.reaction];
					}
				});
			}
		}
	});