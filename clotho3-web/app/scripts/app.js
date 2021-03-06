'use strict';

//full default angular stack
angular.module('clotho.fullPackage', [
    'clotho.foundation',    //core components as described above
    'clotho.commandbar',        //command bar
    'clotho.webapp',        //general web app components
    //additional webapp modules
    'clotho.editor',
    'clotho.interface',
    'clotho.trails',
    'clotho.dna',         //in here temporarily until server can handle all this
    'clotho.construction'
]);

//web application set up
angular.module('clothoRoot', ['clotho.fullPackage'])
.config(function ($routeProvider, $locationProvider) {

    /*
    simulate legacy browser not supporting pushstate (add $provide to DI clause above)
    $provide.decorator('$sniffer', function($delegate) {
        $delegate.history = false;
        return $delegate;
    });
    */

    $locationProvider
        .html5Mode(false)
        .hashPrefix('!');

    $routeProvider
    .when('/', {
        templateUrl: 'views/home.html',
        controller: 'HomeCtrl',
        title : 'Home'
    })
    .when('/login', {
        templateUrl: 'views/login.html',
        controller: 'loginCtrl' //from command module
    })
    .when('/settings', {
        templateUrl: 'views/settings.html',
      controller: 'SettingsCtrl',
        title : 'Settings'
    })
    .when('/about', {
      templateUrl: 'views/about.html',
        title : 'About'
    })
    .when('/team', {
      templateUrl: 'views/team.html',
      controller: 'TeamCtrl',
        title : 'Team'
    })
    .when('/browser', {
      templateUrl: 'views/browser.html',
      controller: 'BrowserCtrl',
        title : 'Browser'
    })
    .when('/edit', {
        redirectTo: '/editor'
    })
    /*
    .when('/editor/query/:queryTerm', {
        templateUrl: 'views/editor.html',
        controller: 'EditorCtrl',
        resolve : {
            queryResult : ['Clotho', '$route', '$q', function (Clotho, $route, $q) {
                var query = '';
                try {
                    query = angular.fromJson($route.current.params.queryTerm);
                } catch (err) {
                    console.warn('editor query malformed: ' + err)
                }

                console.log('querying for editor: ', query);

                if (!angular.isEmpty(query)) {
                    return Clotho.query(query).then(function (data) {
                        return data;
                    });
                } else {
                    return $q.when();
                }
            }],
            deps : ['codemirrorLoader', function(loader) {
                return loader.loadMain();
            }]
        }
    })
    .when('/editor/query', {
        redirectTo : '/editor'
    })
    */
    .when('/editor', {
      templateUrl: 'views/editor.html',
      controller: 'EditorCtrl',
        title : 'Editor',
        reloadOnSearch: false,
        resolve : {
            deps : ['codemirrorLoader', function(loader) {
                return loader.loadMain();
            }]
        }
    })
    .when('/executor', {
        templateUrl: 'views/executor.html',
        controller: 'ExecutorCtrl',
        title : 'Function Executor',
        reloadOnSearch: false
    })
    .when('/trails', {
      templateUrl: 'views/trails.html',
      controller: 'TrailsCtrl',
        title : 'Trails'
    })
    .when('/trail', {
        templateUrl: 'views/trail.html',
        controller: 'TrailCtrl',
        title : 'Trail',
        reloadOnSearch: false,
        resolve : {
            trail : ['Clotho', '$q', '$http', '$route', '$location', 'Trails', function (Clotho, $q, $http, $route, $location, Trails) {
                if (angular.isUndefined($route.current.params.id)) {
                    return $q.when(null);
                }
                var deferred = $q.defer();
                Clotho.get($route.current.params.id).then(function(result) {
                    Trails.compile(result).then(function (compiled) {
                        $route.current.$$route.title = result.name;
                        deferred.resolve(compiled);
                    });
                }, function () {
                    $location.path('/trails')
                });
                return deferred.promise;
            }]
        },
        hotkeys : [
            ['alt+left', 'Previous page of Trail', 'prev()'],
            ['alt+right', 'Next page of Trail', 'next()']
        ]
    })
    .when('/trail-splash', {
        templateUrl: 'views/trail-splash.html',
        controller: 'TrailSplashCtrl',
        title : 'Trail Splash'
    })
    .when('/import', {
        templateUrl:'views/import/intro.html',
        controller : 'ImportCtrl',
        title : 'Import Wizard'
    })
    .when('/import/youtubePlaylist', {
        templateUrl: 'views/import/youtubePlaylist.html',
        controller: 'YoutubePlaylistImportCtrl'
    })
    .when('/import/ape', {
        templateUrl: 'views/import/ape.html',
        controller: 'ApeImportCtrl',
        title: 'GenBank Import'
    })
    .when('/import/ncbi', {
        templateUrl: 'views/import/ncbi.html',
        controller: 'NCBIImportCtrl',
        title: 'NCBI Import'
    })
    .when('/import/facebook', {
        templateUrl: 'views/import/facebook.html',
        controller: 'FacebookImportCtrl',
        title: 'Facebook Import'
    })
  .when('/import/csv', {
    templateUrl: 'views/import/csv.html',
    controller: 'ImportCSVCtrl',
    title: 'CSV Import'
  })

    //testing
    .when('/widgets', {
      templateUrl: 'views/widgets.html',
      controller: 'WidgetsCtrl',
    reloadOnSearch : false
    })

.when('/test/schemaview', {
  templateUrl: 'views/test/schemaview.html',
  controller: 'TestSchemaviewCtrl',
    hotkeys : [
        ['m', 'Show Programmatic Modal', 'createModal()']
    ]
})
.when('/test/quiz', {
  templateUrl: 'views/test/quiz.html',
  controller: 'TestQuizCtrl'
})
.when('/test/construction', {
  templateUrl: 'views/test/construction.html',
  controller: 'TestConstructionCtrl',
    title : 'Construction Test'
})
.when('/test/constructionTrail', {
  templateUrl: 'views/test/contstructiontrail.html',
  controller: 'TestContstructiontrailCtrl',
    title : 'Construction Trail Test',
    reloadOnSearch: false,
    resolve : {
        trail : ['$q', '$http', '$route', 'Trails', function ($q, $http, $route, Trails) {
            var deferred = $q.defer();
                //Clotho.get('org.clothocad.trails.constructionFiles')
                $http.get('models/org.clothocad.trails.constructionFiles.json')
                .then(function(data) {
                    Trails.compile(data.data).then(function (compiled) {
                        deferred.resolve(compiled);
                    });
            });
            return deferred.promise;
        }]
    }
})
.when('/test/focus', {
  templateUrl: 'views/test/focus.html',
  controller: 'TestFocusCtrl'
})
    .when('/test/playground', {
        templateUrl: 'views/test/playground.html'
    })
    .otherwise({
        redirectTo:'/'
    });

})
.run(function($rootScope, $route, $window, interfaceConfig) {

    /****** Config *****/

    $rootScope.$on('$routeChangeSuccess', function(event, current, previous) {
        var title = angular.isDefined(current.$$route) ? current.$$route.title : null;
        //can't use interpolation in document title because ng-app is within body
        $window.document.title = 'Clotho' + (angular.isDefined(title) ? ' | ' + title : '');
    });

    $rootScope.interfaceConfig = interfaceConfig;

});

/*
 angular.element(document).ready(function() {
 angular.bootstrap(document, ['clothoRoot']);
 });
 */
