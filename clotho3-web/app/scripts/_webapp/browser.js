'use strict';

angular.module('clotho.webapp').controller('BrowserCtrl',
function($scope, Clotho, $filter) {

	//todo - pending this working... see GH#410
	Clotho.recent().then(function(result) {
		$scope.resultArray = result;
		$scope.sort(false);
	});

	//todo - store on server and pull
	$scope.savedQueries = [
		{
			name : "Schemas",
			query : {
				schema: "org.clothocad.core.schema.Schema"
			}
		},
		{
			name : "Contains pBca",
			query : {
				name: "{$regex : '[pBca]', $options : 'gi'}"
			}
		}
	];

	$scope.newQuery = {};
	$scope.testNewQuery = function () {
		Clotho.query($scope.newQuery).then(function (results) {
			$scope.newQueryResults = results;
		});
	};
	$scope.saveNewQuery = function () {
		//todo - update once saving user queries on server
		$scope.savedQueries.push($scope.newQuery);
		$scope.newQuery = {};
		$scope.newQueryResults = '';
	};


	//todo - update pending collection schema update - GH#411
	$scope.demoCollection = {
		name : "My Collection",
		author : "uniqueUserID",
		description : "my Collection of all the things.",
		items : [
			{
				id : "myReferentObjectId",
				note : "Note relevant to this collection"
			}
		]
	};

	$scope.setCurrentQuery = function(value) {
		$scope.currentQuery = value;
		Clotho.query(value).then(function (result) {
			$scope.resultArray = result;
		});
	};

	$scope.sort = function(byCat) {
		if (byCat) {
			if (!!$scope.catSort) {return;}

			$scope.catSort = true;
			$scope.recent = $filter('categorize')($scope.recent_array, 'type');
			$scope.recent['Instance'] = $filter('categorize')($scope.recent['Instance'], 'schema.name');
		} else {
			if (!$scope.catSort && angular.isDefined($scope.catSort)) {return;}

			$scope.catSort = false;
			$scope.recent = {
				"entries" : angular.copy($scope.recent_array)
			}
		}
	};
});

	[
		{
			name : "Schemas",
			query : {
				schema: "org.clothocad.core.schema.Schema"
			}
		},
		{
			name : "Contains pBca",
			query : {
				name: "{$regex : '[pBca]', $options : 'gi'}"
			}
		}
	]