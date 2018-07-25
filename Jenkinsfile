#!groovy​

// The 'xchem' Fragalysis Backend Jenkinsfile.

pipeline {

  agent { label 'buildah-slave' }

  environment {
    // Registry details
    USER = 'jenkins'
    REGISTRY = 'docker-registry.default:5000'
    STREAM_IMAGE = "${REGISTRY}/fragalysis-cicd/fragalysis-backend:latest"

    // Slack channel to be used for errors/failures
    SLACK_ALERT_CHANNEL = 'dls-alerts'
  }

  stages {

    stage('Inspect') {
      steps {
          echo "Inspecting..."
      }
    }

    stage('Build Image') {
      steps {
        echo "Building fragalysis-backend..."
        sh "buildah bud --format docker -f Dockerfile -t ${STREAM_IMAGE} ."
      }
    }

    stage('Push Image') {
      steps {
        script {
          TOKEN = sh(script: 'oc whoami -t', returnStdout: true).trim()
        }
        sh "podman login --tls-verify=false --username ${env.USER} --password ${TOKEN} ${env.REGISTRY}"
        sh "buildah push --tls-verify=false ${env.STREAM_IMAGE} docker://${env.STREAM_IMAGE}"
        sh "podman logout ${env.REGISTRY}"
      }
    }

  }

  // Post-job actions.
  // See https://jenkins.io/doc/book/pipeline/syntax/#post
  post {

    failure {
      slackSend channel: "#${SLACK_ALERT_CHANNEL}",
              color: 'danger',
              message: "Fragalysis-Backend build ${env.BUILD_NUMBER} - failed (${env.BUILD_URL})"
    }

    fixed {
      slackSend channel: "#${env.SLACK_ALERT_CHANNEL}",
              color: 'good',
              message: "Fragalysis-Backend build - fixed"
    }

  }

}
