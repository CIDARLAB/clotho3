angular.module('clotho.construction')
/**
 * @ngdoc directive
 * @name constructionFileView
 *
 * @restrict A
 *
 * @description
 * Visualizer for construction files
 */
	.directive('constructionFile', function () {
		return {
			restrict : 'A',
			templateUrl : 'views/_construction/constructionFile.html',
			scope : {
				file : '=constructionFile'
			},
			controller : function constructionFileController ($scope, $element, $attrs) {

				this.getStep = function (index) {
					if (angular.isEmpty($scope.file)) {
						return null;
					}
					return $scope.file.steps[index];
				};

			},
			link : function (scope, element, attrs) {

			}
		}
	});