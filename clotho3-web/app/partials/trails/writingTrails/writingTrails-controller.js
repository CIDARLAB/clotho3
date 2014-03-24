'use strict';

$clotho.extensions.controller('writingTrailsCtrl', function($compile, $scope, Clotho, $http, $location) {
    $scope.getSchema = function(schemaName) {
      Clotho.query({schema: "ClothoSchema", name: schemaName})
        .then(function(results) {
          $scope.demoSchema = results[0]
        });
    };

    $scope.getTrail = function(trailName) {
      Clotho.query({schema: 'Trail', name: trailName})
        .then(function(results) {
          $scope.demoTrail = results[0]
        });
    };

    //temporary hack while Trail schema can't be imported
    $scope.getTrailHack = function() {
        $http.get('partials/trails/writingTrails/demo_trailSchema.json').then(function(data) {
            $scope.demoSchema = data.data;
        });
    };

	//another dumb hack (Because its easier this way)
	$scope.PageExample = {
		title : "Title of Page",
		icon : "exercise",
		controller : "NameOfCustomAngularController",
		mixin : "URL/of/custom/angular/component",
		script: "script to run before compile page",
		onload: "script to run after compile page",
		css : "CSS file to load for this page",
		dictionary : {
			object : "defining keys for interpolation on the page",
			which: "will extend the scope"
		},
		contents : [
			{
				type : "text",
				params : "content will appear on the page in the order it is defined in the contents array. Text is either simple text, or HTML with HTML markup and angular bindings in the form <em>&lt;span ng-non-bindable&gt;{{ variable }}&lt;/span&gt;</em>"
			},
			{
				"type" : "template",
				"params" : "url/to/template/to/include.html"
			},
			{
				type : "wiki",
				params : "A bunch of MediaWiki format text, which will be parsed to HTML"
			},
			{
				type : "markdown",
				params : "A bunch of Markdown format text, which will be parsed to HTML by Showdown.js"
			},
			{
				"type" : "video",
				"params" : {
					id : "10digitYT_ID",
					mini: true,
					autoplay: false,
					params: {
						"height" : 394,
						"width" : 700,
						"ALTERNATIVELY" : "just define params as a string of the videoId or URL to the video, and we'll use all the defaults"
					}

				}
			},
			{
				"type" : "quiz",
				"params" : {
					"title" : "Quiz Title",
					"type" : "See quiz schema for available types",
					"question" : "quiz question, see quiz schema",
					"dictionary" : {},
					"retry" : {}
				}
			},
			{
				type : "hint",
				params : "You can put a little hint button on the right side of the page with a popover displaying the text you put here"
			},
		]
	};

	$scope.$watch(function() {
			return $location.path()
		}, function() {
		$scope.demoSchema = {};
		$scope.demoTrail = {};
	})

})
.controller('writingTrailsUICtrl', function(Clotho, $scope, $focus, $modal) {
    $scope.videoDemo = {
        id : '<11 char videoID>',
        params : '<youtubeParamsObj or URL>',
        autoplay : 'true|false',
        mini : 'true|false'
    };

    $scope.demoSize = {
        "params" : {
            "height" : 394,
            "width" : 700
        }
    };

    $scope.showSearchTypeOut = function () {
        $focus.typeOutSearch('Demo search text');
    };

    $scope.showDialog = function() {
        $modal.messageBox('Demo Modal', 'This is some simple text. It can contain <code>HTML</code>', [{label: "OK", cssClass: "btn-primary", result: true}])
	        .result
          .then(function() {
	          return $modal.messageBox('Chaining', 'The $modal and $focus services return promises so you can chain together UI elements based on their completion / response', [{label: "OK", cssClass: "btn-primary", result: true}])
          })
	        .result
          .then(function() {
              return $focus.typeOutSearch("$focus demo")
          })

    };

    //demo video modal
    $scope.openVideo = function() {
        $modal.video('ivif214cQLE', {});
    };

    $scope.schemas = [];
    Clotho.query({"schema": "Schema"}).then(function(data) {
        $scope.schemas = data;
    });


    $scope.dropdownItems = [
        "The first choice!",
        "And another choice for you.",
        "but wait! A third!"
    ];


})
.controller('writingTrailsClothoCtrl', function(Clotho, $scope) {

});