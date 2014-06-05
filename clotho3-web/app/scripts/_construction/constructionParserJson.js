angular.module('clotho.construction')
/**
 * @ngdoc directive
 * @name constructionParserJson
 *
 * @restrict A
 *
 * @description
 * GUI which takes wet lab construction file and creates Clotho Construction file
 *
 * @example
 *
 */
	.directive('constructionParserJson', function (Clotho) {
		return {
			restrict : 'A',
			scope: {
				wetlab : '=constructionParserJson'
			},
			templateUrl : 'views/_construction/constructionParserJson.html',
			link : function (scope, element, attrs) {
				scope.$watch('wetlab', function (newval) {
					console.log(newval);
					if (newval && newval.length) {

						Clotho.run('clotho.functions.construction.ParseWetlabConstructionFile', [newval])
						.then(function (parsed) {
							console.log(parsed);
							scope.json = parsed;
						});
					}
				})
			}
		}
	});